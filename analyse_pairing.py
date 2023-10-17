import tfs
from matplotlib import pyplot as plt
import numpy as np

BINS = 40

STEP = 0.1
BLUE_INNER = "#5080F0"
BLUE_OUTER = "#000080"

print("=== Sum ================================")
summ = tfs.read("summ_sum.tfs")

print("CORR")
print(summ["CORR_AFTER"]/summ["CORR"])
print("BBEAT")
print(summ["BBEAT_AFTER"]/summ["BBEAT"])

data = summ["CORR_AFTER"]/summ["CORR"]

data_over = [x for x in data if x > 1]
print("CORR > 0")
print(len(data_over)/len(data))

data_under = [x for x in data if x < 1]
print("CORR < 0")
print(len(data_under)/len(data)) 

plt.hist(data_under,
         label="corr (improving)",
         color=BLUE_INNER,
         edgecolor=BLUE_OUTER,
         histtype="stepfilled",
         bins=np.arange(0, 1.1, STEP))

plt.hist(data_over,
         label="corr (deteriorating)",
         color="#F05050",
         edgecolor="#800000",
         histtype="stepfilled",
         bins=np.arange(1, 8, STEP))
plt.xlabel("relative $\\beta$ beating")
plt.ylabel("Frequency")
plt.legend()
plt.savefig("corr.pdf")

plt.hist(summ["BBEAT_AFTER"]/summ["BBEAT"], label="BBEAT", bins=BINS)
plt.legend()
plt.savefig("bbeat.pdf")

print(f"improv: {len(data_under)}")
print(f"deteri: {len(data_over)}")

print("=== Dif ================================")
summ = tfs.read("summ_diff.tfs")

data = summ["CORR_AFTER"]/summ["CORR"]

plt.close()

data_over = [x for x in data if x > 1]
print("CORR > 0")
print(len(data_over)/len(data))

data_under = [x for x in data if x < 1]
print("CORR < 0")
print(len(data_under)/len(data)) 

plt.hist(data_under,
         label="CORR",
         color="#5050F0",
         edgecolor="#000080",
         histtype="stepfilled",
         bins=np.arange(0, 1.1, STEP))

plt.hist(data_over,
         label="CORR",
         color="#F05050",
         edgecolor="#800000",
         histtype="stepfilled",
         bins=np.arange(1, 8, STEP))
plt.legend()
plt.savefig("diff.pdf")

print(summ["CORR_AFTER"]/summ["CORR"])
print(summ["BBEAT_AFTER"]/summ["BBEAT"])

print(np.arange(0, 1.1, 0.2))
