import matplotlib.pyplot as plt
import csv
first = 1130.0
end = 1380.0
diff = end - first

with open('em.csv', newline='') as csvfile:
    r = csv.reader(csvfile, delimiter=';', quotechar='|')
    x = []
    y = []
    counter = 0
    for row in r:
        if counter % 1 == 0 and float(row[0]) > first and float(row[0]) < end and (counter % 100 == 0):
            x.append(float(row[0]))
            y.append(float(row[2]))
        counter += 1
    print(counter)

plt.plot(x, y)
plt.show()

#with open('em.csv', newline='') as csvfile:
#    with open('new.csv', 'w') as new:
#        r = csv.reader(csvfile, delimiter=';', quotechar='|')
#        for row in r:
#            if float(row[0]) < first:
#                counter = 0
#                for k in row:
#                    new.write(k)
#                    if counter != len(row) - 1:
#                        new.write(";")
#                    counter += 1
#                new.write("\n")
#            elif float(row[0]) >= first and float(row[0]) <= end:
#                pass
#            else:
#                pass
#                row[0] = str(float(row[0]) - float(diff))
#                counter = 0
#                for k in row:
#                    new.write(k);
#                    if counter != len(row) - 1:
#                        new.write(";")
#                    counter += 1
#                new.write("\n")

