import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class RNAfolding {

    static class BasePair {
        char firstBase;
        char secondBase;
        int firstIndex;
        int secondIndex;

        BasePair(char firstBase, char secondBase, int firstIndex, int secondIndex) {
            this.firstBase = firstBase;
            this.secondBase = secondBase;
            this.firstIndex = firstIndex;
            this.secondIndex = secondIndex;
        }
    }

    public static void main(String[] cmdArgs) throws IOException {
        if (cmdArgs.length != 1) {
            System.out.println("Usage: RNAFolding <input file name>");
            return;
        }

        String rnaSequence = "";
        int[][] dynamicProgrammingArray = null;
        List<BasePair> pairsList = null;
        int count = 1;
        List<String> sequencesList = new ArrayList<>();

        // Read sequences from the input file
        try (BufferedReader fileReader = new BufferedReader(new FileReader(cmdArgs[0]))) {
            String line;
            while ((line = fileReader.readLine()) != null) {
                if (line.startsWith("**RNA")) {
                    sequencesList.add(fileReader.readLine());
                }
            }
        }

        for (String sequence : sequencesList) {
            rnaSequence = sequence;
            int sequenceLength = rnaSequence.length();
            dynamicProgrammingArray = new int[sequenceLength + 1][sequenceLength + 1];
            pairsList = new ArrayList<>();

            // Perform dynamic programming to find the optimal secondary structure
            for (int length = 5; length < sequenceLength; length++) {
                for (int start = 1; start <= sequenceLength - length; start++) {
                    int end = start + length;
                    dynamicProgrammingArray[start][end] = dynamicProgrammingArray[start][end - 1];

                    for (int mid = start; mid < end - 4; mid++) {
                        char leftChar = rnaSequence.charAt(mid - 1);
                        char rightChar = rnaSequence.charAt(end - 1);

                        if (areComplementary(leftChar, rightChar)) {
                            dynamicProgrammingArray[start][end] = Math.max(
                                    dynamicProgrammingArray[start][end],
                                    dynamicProgrammingArray[start][mid - 1] + dynamicProgrammingArray[mid + 1][end - 1]
                                            + 1);
                        }
                    }
                }
            }

            // Output the optimal secondary structure and base pairs
            System.out
                    .println("** RNA-" + count + ", length=" + rnaSequence.length() + ", Optimal secondary structure:");
            findOptimalBasePairs(1, sequenceLength, dynamicProgrammingArray, rnaSequence, pairsList);

            pairsList.sort((pair1, pair2) -> {
                if (pair1.firstIndex != pair2.firstIndex) {
                    return Integer.compare(pair1.firstIndex, pair2.firstIndex);
                } else {
                    return Integer.compare(pair1.secondIndex, pair2.secondIndex);
                }
            });

            for (BasePair pair : pairsList) {
                System.out.println(
                        pair.firstBase + "-" + pair.secondBase + " (" + pair.firstIndex + "," + pair.secondIndex + ")");
            }

            System.out.println("Total number of base pairs: " + dynamicProgrammingArray[1][sequenceLength]);
            System.out.println();
            count++;
        }

        System.out.println("*** Asg 5 by Drithi Madagani***");
    }

    private static boolean areComplementary(char a, char b) {
        // Check if two RNA bases are complementary (e.g., A-U or C-G)
        return (a == 'A' && b == 'U') || (a == 'U' && b == 'A') || (a == 'C' && b == 'G') || (a == 'G' && b == 'C');
    }

    private static void findOptimalBasePairs(int i, int j, int[][] dp, String sequence, List<BasePair> pairs) {
        if (i < j) {
            if (dp[i][j] == dp[i][j - 1]) {
                findOptimalBasePairs(i, j - 1, dp, sequence, pairs);
            } else {
                int maxBasePairs = dp[i][j];
                int splitIndex = -1;

                for (int k = i; k < j - 4; k++) {
                    if (areComplementary(sequence.charAt(k - 1), sequence.charAt(j - 1))
                            && dp[i][j] == dp[i][k - 1] + dp[k + 1][j - 1] + 1) {
                        splitIndex = k;
                        break;
                    }
                }

                if (splitIndex != -1) {
                    pairs.add(new BasePair(sequence.charAt(splitIndex - 1), sequence.charAt(j - 1), splitIndex, j));
                    findOptimalBasePairs(i, splitIndex - 1, dp, sequence, pairs);
                    findOptimalBasePairs(splitIndex + 1, j - 1, dp, sequence, pairs);
                }
            }
        }
    }
}
