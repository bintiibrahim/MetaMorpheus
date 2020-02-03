using EngineLayer;
using EngineLayer.FdrAnalysis;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace TaskLayer
{
    public static class ProteinProbability
    {
        public static PostSearchAnalysisParameters Parameters { get; set; }
        public static CommonParameters CommonParameters { get; set; }

        public static void CalculateProteinProbabilities(PostSearchAnalysisParameters parameters, CommonParameters commonParameters, List<PeptideSpectralMatch> peptides)
        {
            Parameters = parameters;
            CommonParameters = commonParameters;

            // commented out for now because this step is done in post search analysis
            // 1. filter to get best (highest-scoring) PSM per peptide (copied from post-search analysis)
            //List<PeptideSpectralMatch> peptides = Parameters.AllPsms.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList(); // mod peptides are treated as different peptides

            //new FdrAnalysisEngine(peptides, Parameters.NumNotches, CommonParameters, new List<string> { Parameters.SearchTaskId }, "Peptide").Run();

            //if (!Parameters.SearchParameters.WriteDecoys)
            //{
            //    peptides.RemoveAll(b => b.IsDecoy);
            //}
            //if (!Parameters.SearchParameters.WriteContaminants)
            //{
            //    peptides.RemoveAll(b => b.IsContaminant);
            //}
            //peptides.RemoveAll(p => p.FdrInfo.QValue > CommonParameters.QValueOutputFilter);

            // 2. populate dictionary for protein-peptide associations
            var peptidesForProteins = new Dictionary<Protein, List<PeptideSpectralMatch>>();
            var proteinsForPeptides = new Dictionary<string, List<Protein>>();

            foreach (PeptideSpectralMatch bestScoringPsm in peptides)
            {
                foreach (Protein protein in bestScoringPsm.BestMatchingPeptides.Select(p => p.Peptide.Protein))
                {
                    if (peptidesForProteins.TryGetValue(protein, out List<PeptideSpectralMatch> peptidesForThisProtein))
                    {
                        peptidesForThisProtein.Add(bestScoringPsm);
                    }
                    else
                    {
                        peptidesForProteins.Add(protein, new List<PeptideSpectralMatch> { bestScoringPsm });
                    }


                    if (proteinsForPeptides.TryGetValue(GetSequence(bestScoringPsm), out List<Protein> proteinsForThisPeptide))
                    {
                        proteinsForThisPeptide.Add(protein);
                    }
                    else
                    {
                        proteinsForPeptides.Add((GetSequence(bestScoringPsm)), new List<Protein> { protein });
                    }
                }
            }

            //// calculate NSP for all peptides
            //foreach (var kvp in peptidesForProteins.Where(p => p.Value.Count > 1))
            //{
            //    var peptidesForThisProtein = new List<PeptideSpectralMatch>(kvp.Value);
            //    foreach (var peptide in kvp.Value)
            //    {
            //        var neighborPeptides = new List<PeptideSpectralMatch>(kvp.Value);
            //        neighborPeptides.Remove(peptide);

            //        peptide.NSP = neighborPeptides.Select(n => n.FdrInfo.PEP).Aggregate(0.0, (acc, val) => acc + val);
            //    }
            //}

            // 3. compute protein probabilities for proteins with unique peptides only (for now)
            var proteinProbabilities = new Dictionary<Protein, double>();
            foreach (var kvp in peptidesForProteins)
            {
                // 3a. exclude shared peptides
                var uniquePeptidesForThisProtein = kvp.Value.Where(p => proteinsForPeptides[GetSequence(p)].Count == 1);

                // 3b. find probability of incorrect peptides for this protein i.e. the product of PEPs
                double probabilityIncorrectPeptides = uniquePeptidesForThisProtein.Select(p => p.FdrInfo.PEP).Aggregate(1.0, (acc, val) => acc * val);

                // 3c. find the probability of protein i.e. at least 1 correct peptide
                double proteinProbabiity = 1 - probabilityIncorrectPeptides;

                proteinProbabilities.Add(kvp.Key, proteinProbabiity);
            }

            // 4. Write results to file
            WriteToFile(Path.Combine(Parameters.OutputFolder, "ProteinProbabilities"), proteinProbabilities, peptidesForProteins);

        }

        private static void WriteToFile(string filePath, Dictionary<Protein, double> proteinProbabilites, Dictionary<Protein, List<PeptideSpectralMatch>> peptidesForProteins)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                foreach (var kvp in proteinProbabilites)
                {
                    string peptidesForThisProtein = string.Join(",", peptidesForProteins[kvp.Key].Select(p => GetSequence(p)));
                    output.WriteLine(string.Join("\t", new string[] { kvp.Key.ToString(), peptidesForThisProtein, kvp.Value.ToString() }));
                }
            }
        }

        private static string GetSequence(PeptideSpectralMatch psm)
        {
            return psm.FullSequence ?? string.Join("|", psm.BestMatchingPeptides.Select(pep => pep.Peptide.FullSequence));
        }

    }
}
