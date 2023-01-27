import csv
import math
from scipy import stats

residuals = []
times = []

first_row = True
with open('residual_1.csv') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        if first_row:
            first_row = False
        else:
            tmp = [float(row[2]), 0.0, 0.0, 0.0]
            residuals.append(tmp)
            times.append(float(row[1]))

first_row = True
i = 0
with open('residual_2.csv') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        if first_row:
            first_row = False
        else:
            residuals[i][1] = float(row[2])
            i = i + 1

first_row = True
i = 0
with open('residual_3.csv') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        if first_row:
            first_row = False
        else:
            residuals[i][2] = float(row[2])
            i = i + 1

first_row = True
i = 0
with open('residual_4.csv') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        if first_row:
            first_row = False
        else:
            residuals[i][3] = float(row[2])
            i = i + 1

h_logs = [math.log(1.0 / 3.0), math.log(2.0 / 7.0), math.log(1.0 / 4.0), math.log(2.0 / 9.0)]
h_logs = [math.log(0.110126), math.log(0.0957201), math.log(0.0826757), math.log(0.0755121)]
orders = []
first_row = True
i = 0
for res in residuals:
    if first_row:
        first_row = False
    else:
        res_logs = [math.log(res[0]), math.log(res[1]), math.log(res[2]), math.log(res[3])]
        fit = stats.linregress(h_logs, res_logs)
        # fit = stats.linregress(res_logs, h_logs)
        orders.append([times[i], fit.slope])
    i = i + 1

with open('orders.csv', 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Time", "Order"])
    for o in orders:
        csvwriter.writerow(o)
