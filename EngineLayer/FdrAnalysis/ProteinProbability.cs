using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    public static class ProteinProbability
    {
        // for testing purposes i have this method called from FdrAnalysisEngine
        // but only works if allPsms count > 100 to use PEP
        public static void CalculateProteinProbabilities(List<PeptideSpectralMatch> allPsms)
        {
            var peptidesToProtein = new Dictionary<Protein, List<string>>();
            var proteinsToPeptide = new Dictionary<string, List<Protein>>();
            var proteinProbabiltiies = new Dictionary<Protein, double>();

            foreach (var (Notch, Peptide) in allPsms.SelectMany(p => p.BestMatchingPeptides).Distinct())
            {
                string sequence = Peptide.BaseSequence;

                // 1. create dictionary with protein as keys that correspond to a list of peptide sequences
                if (peptidesToProtein.TryGetValue(Peptide.Protein, out List<string> peptides))
                {
                    peptides.Add(sequence);
                }
                else
                {
                    peptidesToProtein.Add(Peptide.Protein, new List<string> { sequence });
                }


                // 2. create dictionary with peptide sequence as keys that correspond to a list of proteins
                if (proteinsToPeptide.TryGetValue(sequence, out List<Protein> proteins))
                {
                    proteins.Add(Peptide.Protein);
                }
                else
                {
                    proteinsToPeptide.Add(sequence, new List<Protein> { Peptide.Protein });
                }
            }

            var PEPToPeptide = new Dictionary<string, List<double>>();
            foreach (var kvp in peptidesToProtein)
            {
                foreach (var sequence in kvp.Value)
                {
                    // 3. find psms that contain unique peptides associated with this protein to access PEP values from psm.FdrInfo.PeptidePEP
                    var associatedPsms = allPsms.Where(psm => psm.BestMatchingPeptides.Select(p => p.Peptide)
                        .Any(p => p.BaseSequence.Equals(sequence) 
                        && p.Protein == kvp.Key
                        && proteinsToPeptide[sequence].Count == 1));

                    foreach (var psm in associatedPsms)
                    {
                        var associatedPeptides = psm.BestMatchingPeptides.Select(p => p.Peptide)
                            .Where(p => p.Protein == kvp.Key
                            && p.BaseSequence.Equals(sequence)
                            && proteinsToPeptide[sequence].Count == 1);

                        foreach (var peptide in associatedPeptides)
                        {
                            var PEP = psm.FdrInfo.PeptidePEP[peptide];

                            if (PEPToPeptide.TryGetValue(sequence, out List<double> PEPs))
                            {
                                PEPs.Add(PEP);
                            }
                            else
                            {
                                PEPToPeptide.Add(sequence, new List<double>() { PEP });
                            }
                        }
                    }
                }

                // 4. compute protein probability as the probability of at least 1 correct peptide assignment
                double probability = 1;
                foreach (var pep_kvp in PEPToPeptide)
                {
                    probability *= pep_kvp.Value.Max();
                }

                double proteinProbability = 1 - probability;
                proteinProbabiltiies.Add(kvp.Key, proteinProbability);
            }
        }
    }
}
