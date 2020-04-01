using Proteomics;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer
{
    public class ProteinData
    {
        public static readonly string[] featuresForTraining = new string[] {
            "AveragePEP", "PercentageOfUniquePeptides", "TotalPeptideCount", "PercentageSequenceCoverage",
            "NumProteinAccessions", "BestPeptideScore", "LongestPepSeries" };

        public string ProteinGroupName { get; set; }

        public float AveragePEP { get; set; }

        public float PercentageOfUniquePeptides { get; set; }

        public float TotalPeptideCount { get; set; }

        public float PercentageSequenceCoverage { get; set; }

        public float NumProteinAccessions { get; set; }

        public float BestPeptideScore { get; set; }

        public float LongestPepSeries { get; set; }

        public bool Label { get; set; }

        public bool IsDecoy { get; set; }

        public bool IsContaminant { get; set; }

        public double QValue { get; set; }

        public string GetTabSeparatedHeader()
        {
            var sb = new StringBuilder();
            sb.Append("Protein Accession" + '\t');
            sb.Append("Number of Proteins" + '\t');
            sb.Append("Total Peptide Count" + '\t');
            sb.Append("Percentage of Unique Peptides" + '\t');
            sb.Append("Normalized Longest Peptide Series" + '\t');
            sb.Append("Sequence Coverage Fraction" + '\t');
            sb.Append("Best Peptide Score" + '\t');
            sb.Append("Average PEP" + '\t');
            sb.Append("Q-Value" + '\t');
            sb.Append("Protein Decoy/Contaminant/Target" + '\t');
            sb.Append("Label" + '\t');
            return sb.ToString();
        }

        public override string ToString()
        {
            var sb = new StringBuilder();

            // list of protein accessions
            sb.Append(ProteinGroupName);
            sb.Append("\t");

            // number of proteins in group
            sb.Append(NumProteinAccessions);
            sb.Append("\t");

            // total peptides
            sb.Append(TotalPeptideCount);
            sb.Append("\t");

            // percentage of unique peptides
            sb.Append(PercentageOfUniquePeptides);
            sb.Append("\t");

            // longest peptide series, normalized
            sb.Append(LongestPepSeries);
            sb.Append("\t");

            // sequence coverage fraction
            sb.Append(PercentageSequenceCoverage);
            sb.Append("\t");

            // best score among all proteins in group
            sb.Append(BestPeptideScore);
            sb.Append("\t");

            // average PEP value in group
            sb.Append(AveragePEP);
            sb.Append("\t");

            // q-value
            sb.Append(QValue);
            sb.Append("\t");

            // decoy/contaminant/target
            if (IsDecoy)
            {
                sb.Append("D");
            }
            else if (IsContaminant)
            {
                sb.Append("C");
            }
            else
            {
                sb.Append("T");
            }
            sb.Append("\t");

            // training label
            sb.Append(Label);
            sb.Append("\t");

            return sb.ToString();
        }
    }
}
