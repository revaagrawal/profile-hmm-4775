import csv
import os

def main():
    expected_positives = 0
    true_positives = 0
    false_positives = 0

    expected_negatives = 0
    true_negatives = 0
    false_negatives = 0

    file_name = "final_results.csv"

    with open(file_name, "r") as file:
        reader = csv.reader(file)
        for row in reader:
            if len(row) == 0:
                continue
            expected = row[-2]
            match = row[-1]
            print(expected, match)
            if row[-1] != "ERROR":
                if expected == "1":
                    expected_positives += 1
                    if match == "1":
                        true_positives += 1
                    else:
                        false_negatives += 1
                else:
                    expected_negatives += 1
                    if match == "1":
                        true_negatives += 1
                    else:
                        false_positives += 1

    print(expected_positives)
    print(true_positives)
    print(false_positives)

    print(expected_negatives)
    print(true_negatives)
    print(false_negatives)

    print("true positive rate:", true_positives / (true_positives + false_negatives))
    print("true negative rate:", true_negatives / (true_negatives + false_positives))
    print("false positive rate:", false_positives / (false_positives + true_negatives))
    print("false negative rate:", false_negatives / (false_negatives + true_positives))
 
if __name__ == "__main__":
    main()
