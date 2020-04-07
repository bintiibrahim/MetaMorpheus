using EngineLayer.FdrAnalysis;
using Microsoft.ML;
using Microsoft.ML.Data;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using static Microsoft.ML.DataOperationsCatalog;

namespace EngineLayer
{
    public class ProteinProbabilityAnalysis
    {
        private const double PeptidePValueCutoff = 0.0;
        private static List<string> ProteinIDs { get; set; }

    public static string ComputeProteinProbabilities(List<ProteinGroup> proteinGroups, bool modPeptidesAreDifferent, string filePath, List<string> dbFileList, List<string> rawFileList)
        {
            // create ML context
            MLContext mlContext = new MLContext();
            var proteinData = CreateProteinData(proteinGroups, modPeptidesAreDifferent);
            IDataView dataView = mlContext.Data.LoadFromEnumerable(proteinData.AsEnumerable());

            // write training features
            SetProteinIDs(dbFileList, rawFileList);
            WriteProteinDataToTsv(proteinData, filePath);

            // split data into train and test
            TrainTestData trainTestSplit = mlContext.Data.TrainTestSplit(dataView, testFraction: 0.1, null, 42);
            IDataView trainingData = trainTestSplit.TrainSet;
            IDataView testData = trainTestSplit.TestSet;

            // get training features
            var features = ProteinData.featuresForTraining;

            // create model
            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features");

            var pipeline = mlContext.Transforms.Concatenate("Features", features)
                .Append(mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features"));

            // train model
            var trainedModel = pipeline.Fit(trainingData);

            // create prediction engine
            var predictionEngine = mlContext.Model.CreatePredictionEngine<ProteinData, TruePositivePrediction>(trainedModel);

            // predict protein probability for each protein entry
            foreach (ProteinGroup pg in proteinGroups)
            {
                ProteinData pd = CreateProteinDataEntry(pg, modPeptidesAreDifferent, false); // fixme
                var proteinPrediction = predictionEngine.Predict(pd);
                pg.Probability = proteinPrediction.Probability;
            }

            var predictions = trainedModel.Transform(testData);

            CalibratedBinaryClassificationMetrics metrics;
            try
            {
                metrics = mlContext.BinaryClassification.Evaluate(data: predictions, labelColumnName: "Label", scoreColumnName: "Score");
                return PrintBinaryClassificationMetrics(trainer.ToString(), metrics);
            }
            catch
            {
                return "";
            }
        }

        public static string PrintBinaryClassificationMetrics(string name, CalibratedBinaryClassificationMetrics metrics)
        {
            StringBuilder s = new StringBuilder();
            s.AppendLine("************************************************************");
            s.AppendLine("*       Metrics for Determination of PEP Using Binary Classification      ");
            s.AppendLine("*-----------------------------------------------------------");
            s.AppendLine("*       Accuracy:  " + metrics.Accuracy.ToString());
            s.AppendLine("*       Area Under Curve:  " + metrics.AreaUnderRocCurve.ToString());
            s.AppendLine("*       Area under Precision recall Curve:  " + metrics.AreaUnderPrecisionRecallCurve.ToString());
            s.AppendLine("*       F1Score:  " + metrics.F1Score.ToString());
            s.AppendLine("*       LogLoss:  " + metrics.LogLoss.ToString());
            s.AppendLine("*       LogLossReduction:  " + metrics.LogLossReduction.ToString());
            s.AppendLine("*       PositivePrecision:  " + metrics.PositivePrecision.ToString());
            s.AppendLine("*       PositiveRecall:  " + metrics.PositiveRecall.ToString());
            s.AppendLine("*       NegativePrecision:  " + metrics.NegativePrecision.ToString());
            s.AppendLine("*       NegativeRecall:  " + metrics.NegativeRecall.ToString());
            s.AppendLine("************************************************************");
            return s.ToString();
        }

        // build the protein data set that we will feed into the model
        public static List<ProteinData> CreateProteinData(List<ProteinGroup> proteinGroups, bool treatModPeptidesDifferent)
        {
            List<ProteinData> proteinDataList = new List<ProteinData>();

            // due to parsimony, each protein group has 1 protein
            foreach (ProteinGroup pg in proteinGroups)
            {
                bool label;
                if (pg.IsDecoy || pg.QValue > 0.05)
                {
                    label = false;
                    ProteinData newProteinData = CreateProteinDataEntry(pg, treatModPeptidesDifferent, label);
                    proteinDataList.Add(newProteinData);
                }
                else if (!pg.IsDecoy && pg.QValue <= 0.01)
                {
                    label = true;
                    ProteinData newProteinData = CreateProteinDataEntry(pg, treatModPeptidesDifferent, label);
                    proteinDataList.Add(newProteinData);
                }
            }

            return proteinDataList;
        }

        // each entry contains information about the protein and its features
        public static ProteinData CreateProteinDataEntry(ProteinGroup pg, bool treatModPeptidesDifferent, bool label)
        {
            double averagePEP = 0.0;
            double percentageUnique = 0;
            int totalPeptideCount = 0;
            double sequenceCoverageFraction = 0.0;
            int numProteinAccessions = 0;
            double bestScore = 0;

            percentageUnique = (double)pg.UniquePeptides.Count / pg.AllPeptides.Count;
            totalPeptideCount = pg.AllPeptides.Count;
            sequenceCoverageFraction = pg.SequenceCoverageFraction.First();
            numProteinAccessions = pg.Proteins.Count;
            bestScore = pg.BestPeptideScore; 

            List<PeptideSpectralMatch> peptides;
            if (treatModPeptidesDifferent)
            {
                peptides = pg.AllPsmsBelowOnePercentFDR.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
            }
            else
            {
                peptides = pg.AllPsmsBelowOnePercentFDR.GroupBy(b => b.BaseSequence).Select(b => b.FirstOrDefault()).ToList();
            }

            // compute average PEP
            double sumPEP = 0.0;
            foreach (PeptideSpectralMatch peptide in peptides)
            {
                if (peptide.FdrInfo.PEP <= PeptidePValueCutoff) // p value shouldnt be less than 0
                {
                    continue;
                }

                sumPEP += peptide.FdrInfo.PEP;
            }

            averagePEP = sumPEP / peptides.Count;

            // find longest peptide series
            int numPeptidesInSeries = 0;
            int series = 0;

            var oneBasedResiduesByProtein = new Dictionary<Protein, HashSet<(int, int)>>();
            foreach (var psm in pg.AllPsmsBelowOnePercentFDR)
            {
                if (psm.BaseSequence != null)
                {
                    foreach (var pep in psm.BestMatchingPeptides.Select(p => p.Peptide))
                    {
                        if (pg.Proteins.Contains(pep.Protein))
                        {
                            (int, int) uniqueRes = (pep.OneBasedStartResidueInProtein, pep.OneBasedEndResidueInProtein);

                            if (oneBasedResiduesByProtein.TryGetValue(pep.Protein, out HashSet<(int,int)> oneBasedResidues))
                            {
                                oneBasedResidues.Add(uniqueRes);
                            }
                            else
                            {
                                oneBasedResiduesByProtein.Add(pep.Protein, new HashSet<(int, int)> { uniqueRes });
                            }
                        }
                    }
                }
            }

            // weigh unique and shared peptides equally
            var longestPeptideSeries = new List<double>();
            foreach (var prot in pg.Proteins)
            {
                var orderedOneBasedResidues = oneBasedResiduesByProtein[prot].OrderBy(t => t.Item1).ToList(); // order by start residue
                for (int i = 0; i < orderedOneBasedResidues.Count() - 1; ++i)
                {
                    var current = orderedOneBasedResidues[i];
                    var next = orderedOneBasedResidues[i + 1]; // fixme

                    // succeeding peptide match overlaps or comes right after
                    if (next.Item1 - current.Item2 <= 1)
                    {
                        ++numPeptidesInSeries;
                    }
                    else
                    {
                        series = numPeptidesInSeries > series ? numPeptidesInSeries : series; // update longest peptide series
                        numPeptidesInSeries = i + 1 < (orderedOneBasedResidues.Count - 1) ? 0 : numPeptidesInSeries;
                    }
                }

                longestPeptideSeries.Add((double)series / prot.Accession.Length);
            }

            return new ProteinData
            {
                ProteinGroupName = string.Join("|", pg.Proteins.OrderBy(p => p.Accession).Select(p => p.Accession)),
                AveragePEP = (float)averagePEP,
                PercentageOfUniquePeptides = (float)percentageUnique,
                TotalPeptideCount = totalPeptideCount,
                PercentageSequenceCoverage = (float)sequenceCoverageFraction,
                NumProteinAccessions = numProteinAccessions,
                BestPeptideScore = (float)bestScore,
                LongestPepSeries = (float)longestPeptideSeries.Max(),
                IsDecoy = pg.IsDecoy,
                IsContaminant = pg.IsContaminant,
                QValue = pg.QValue,
                Label = label
            };
        }

        private static string GetPrestFileName(string db)
        {
            return Regex.Replace(Path.GetFileNameWithoutExtension(db), @"prest_([a-z0-9]+)_", "");
        }

        private static void SetProteinIDs(List<string> dbFileList, List<string> rawFileList)
        {
            ProteinIDs = new List<string>();

            // get runType
            string runType = Path.GetFileNameWithoutExtension(rawFileList.First()).Replace("mixture", ""); // blank, A, B, or AB
            runType = Regex.Replace(runType, @"rep[1-9]{1}", "");

            var associatedProteinDBs = new List<string>();
            switch (runType)
            {
                case "A":
                    associatedProteinDBs.Add(dbFileList.Where(db => 
                        GetPrestFileName(db).Equals("a")).First());
                    goto case "blank";

                case "B":
                    associatedProteinDBs.Add(dbFileList.Where(db =>
                       GetPrestFileName(db).Equals("b")).First());
                    goto case "blank";

                case "AB":
                    associatedProteinDBs.AddRange(dbFileList.Where(db =>
                       GetPrestFileName(db).Equals("a") ||
                       GetPrestFileName(db).Equals("b")));
                    goto case "blank";

                case "blank":
                    associatedProteinDBs.Add(dbFileList.Where(db =>
                       GetPrestFileName(db).Equals("EColi")).First());
                    break;
            }

            // if runType A, return A & EColi
            // if runType B, return B & EColi
            // if runType AB, return A & B & EColi
            // if runType blank, return EColi
            foreach (var dbFile in associatedProteinDBs)
            {
                string line;
                using (StreamReader file = new StreamReader(dbFile))
                {
                    while ((line = file.ReadLine()) != null)
                    {
                        if (line.StartsWith(">"))
                        {
                            var proteinId = dbFile.Contains("EColi") ? line.Split("|")[1] : line.Substring(1);
                            ProteinIDs.Add(proteinId);
                        }
                    }
                }
            }
        }

        private static void WriteProteinDataToTsv(List<ProteinData> proteinData, string filePath)
        {
            if (proteinData != null && proteinData.Any())
            {
                using (StreamWriter output = new StreamWriter(filePath))
                {
                    output.WriteLine(proteinData.First().GetTabSeparatedHeader());
                    foreach (ProteinData pd in proteinData)
                    {
                        pd.SetKnownLabel(ProteinIDs);
                        output.WriteLine(pd);
                    }
                }
            }
        }

    }
}
