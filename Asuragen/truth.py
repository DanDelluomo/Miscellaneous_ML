from scipy.stats import variation
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

class Assay:
    def assay_analyzer(self, file_to_analyze):

        plt.title('Samples: Scaled Ratios and Copy Numbers')
        plt.xlabel("Sample Number")
        plt.ylabel("Scaled Ratio")
        red_patch = mpatches.Patch(color='red', label='0')
        green_patch = mpatches.Patch(color='green', label='1')
        blue_patch = mpatches.Patch(color='blue', label='2')
        cyan_patch = mpatches.Patch(color='cyan', label='3+')
        plt.legend(title='Expected Copies',handles=[red_patch, green_patch, blue_patch, cyan_patch],loc="upper center")

        plt.axhline(0.364)
        plt.text(1,0.25,'Bin: 0 copies',rotation=0)

        plt.axhline(0.7)
        plt.text(1,0.73,'Bin: 1 copy',rotation=0)
        plt.axhline(0.807)

        plt.axhline(0.814)
        plt.text(1,0.96,'Bin: 2 copies',rotation=0)
        plt.axhline(1.100)

        plt.axhline(1.209)
        plt.text(1,1.33,'Bin: 3+ copies',rotation=0)

        cotnversion_factor = None
        assay_accuracy = []
        zero_copies_ratios = []
        one_copy_ratios = []
        two_copies_ratios = []
        three_plus_copies_ratios = []

        filepath = file_to_analyze
        print("Actual Copy Numbers:")

        with open(filepath) as file:
           line = file.readline()
           count = 0
           while line:
               # ignore the first line of input file
               if count == 0:
                   count += 1
                   line = file.readline()
               else:
                   sample = line.split()
                   if sample[0] == "Calibrator.fsa":
                       conversion_factor = int(sample[3]) / int(sample[2])
                   else:
                       # ignore samples with incomplete data (fewer than 5 data points)
                       if len(sample) > 4:
                           scaled_ratio = int(sample[3]) / int(sample[2]) * (1 / conversion_factor)
                           expected_copies = int(sample[4])

                           if 0 <= scaled_ratio <= 0.364:
                               print(sample[1],": 0")
                               if expected_copies == 0:
                                   assay_accuracy.append(1)
                               else:
                                   assay_accuracy.append(0)
                           elif 0.700 <= scaled_ratio <= 0.807:
                               print(sample[1],": 1")
                               if expected_copies == 1:
                                   assay_accuracy.append(1)
                               else:
                                   assay_accuracy.append(0)
                           elif 0.814 <= scaled_ratio <= 1.100:
                               print(sample[1],": 2")
                               if expected_copies == 2:
                                   assay_accuracy.append(1)
                               else:
                                   assay_accuracy.append(0)
                           elif scaled_ratio >= 1.209:
                               print(sample[1],": 3+")
                               if expected_copies >= 3:
                                   assay_accuracy.append(1)
                               else:
                                   assay_accuracy.append(0)
                           else:
                               print(sample[1],": RERUN")
                               assay_accuracy.append(0)

                           if expected_copies == 0:
                               zero_copies_ratios.append(scaled_ratio)
                               plt.scatter(count, scaled_ratio, color="red")
                           elif expected_copies == 1:
                               one_copy_ratios.append(scaled_ratio)
                               plt.scatter(count, scaled_ratio, color="green")
                           elif expected_copies == 2:
                               two_copies_ratios.append(scaled_ratio)
                               plt.scatter(count, scaled_ratio, color="blue")
                           elif expected_copies >= 3:
                               three_plus_copies_ratios.append(scaled_ratio)
                               plt.scatter(count, scaled_ratio, color="cyan")
                       else:
                          print(sample[1],": RERUN")
                          assay_accuracy.append(0)

                   line = file.readline()
                   count += 1

        print("The coefficient of variation for samples with 0 expected copies is", round(variation(zero_copies_ratios),3))
        print("The coefficient of variation for samples with 1 expected copy is", round(variation(one_copy_ratios),3))
        print("The coefficient of variation for samples with 2 expected copies is", round(variation(two_copies_ratios),3))
        print("The coefficient of variation for samples with 3+ expected copies is", round(variation(three_plus_copies_ratios),3))

        print("The overall accuracy for samples with an expected copy number of 0 is",round(sum(i <= 0.364 for i in zero_copies_ratios) / len(zero_copies_ratios) * 100, 2),"%")
        print("The overall accuracy for samples with an expected copy number of 1 is",round(sum(i >= 0.700 and i <= 0.807 for i in one_copy_ratios) / len(one_copy_ratios) * 100, 2),"%")
        print("The overall accuracy for samples with an expected copy number of 2 is",round(sum(i >= 0.814 and i <= 1.100 for i in two_copies_ratios) / len(two_copies_ratios) * 100, 2),"%")
        print("The overall accuracy for samples with an expected copy number of 3+ is",round(sum(i >= 1.209 for i in three_plus_copies_ratios) / len(three_plus_copies_ratios) * 100, 2),"%")

        print("The overall accuracy of the assay is",round(assay_accuracy.count(1)/len(assay_accuracy) * 100, 2),"%")

        plt.show()

a = Assay()
a.assay_analyzer("inputFile.txt")
