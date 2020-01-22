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

            // 1. filter to get best (highest-scoring) PSM per peptide (copied from post-search analysis task)
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

                    if (proteinsForPeptides.TryGetValue(bestScoringPsm.FullSequence, out List<Protein> proteinsForThisPeptide))
                    {
                        proteinsForThisPeptide.Add(protein);
                    }
                    else
                    {
                        proteinsForPeptides.Add(bestScoringPsm.FullSequence, new List<Protein> { protein });
                    }
                }
            }

            // 4. compute protein probabilities for proteins with unique peptides only (for now)
            var proteinProbabilities = new Dictionary<Protein, double>();
            foreach (var kvp in peptidesForProteins)
            {
                // 4a. exclude shared peptides
                var uniquePeptidesForThisProtein = kvp.Value.Where(peptide => proteinsForPeptides[peptide.FullSequence].Count == 1);

                // 4b. find probability of incorrect peptides for this protein i.e. the product of PEPs
                double probabilityIncorrectPeptides = uniquePeptidesForThisProtein.Select(p => p.FdrInfo.PEP).Aggregate(1.0, (acc, val) => acc * val);

                // 4c. find the probability of protein i.e. at least 1 correct peptide
                double proteinProbabiity = 1 - probabilityIncorrectPeptides;

                proteinProbabilities.Add(kvp.Key, proteinProbabiity);
            }

            // 5. Write to file
            WriteToFile(Path.Combine(Parameters.OutputFolder, "ProteinProbabilities"), proteinProbabilities);

        }

        private static void WriteToFile(string filePath, Dictionary<Protein, double> proteinProbabilites)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine(PeptideSpectralMatch.GetTabSeparatedHeader());
                foreach (var kvp in proteinProbabilites)
                {
                    output.WriteLine(string.Join("\t", new string[] { kvp.Key.ToString(), kvp.Value.ToString() }));
                }
            }
        }

    }
}
