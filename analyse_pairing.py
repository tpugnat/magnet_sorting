import tfs
from matplotlib import pyplot as plt
import numpy as np

BINS = 40

STEP = 0.1
BLUE_INNER = "#5080F0"
BLUE_OUTER = "#000080"

TYPE_SORTING="Sum"
TYPE_SORTING="SumBeta"
TYPE_SORTING="SumBBeating"


print(f"=== {TYPE_SORTING} ================================")
summ = tfs.read(f"summ_{TYPE_SORTING.lower()}.tfs")

data = summ["CORRX_AFTER"]/summ["CORRX"]

data_over = [x for x in data if x >= 1]
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
fig = plt.gcf()
fig.set_size_inches(8, 4)
plt.savefig(f"corr_{TYPE_SORTING.lower()}.pdf")

plt.hist(summ["BBEATX_AFTER"]/summ["BBEATX"], label="BBEATX",
         bins=BINS,
         )
plt.legend()
plt.savefig(f"bbeat_{TYPE_SORTING.lower()}.pdf")

print(f"improv: {len(data_under)}")
print(f"deteri: {len(data_over)}")

plt.close()

# --------------------------------------------------------------------------------------------------

XMAX = 0.03
step = XMAX / BINS

plt.hist(summ["CORRX"],
         label="before sorting",
         color="#F05050",
         edgecolor="#F05050",
         histtype="step",
         linewidth=2,
         bins=np.arange(0, XMAX, step)
         )

plt.hist(summ["CORRX_AFTER"],
         label="after sorting",
         color=BLUE_INNER,
         edgecolor=BLUE_INNER,
         histtype="step",
         linewidth=2,
         bins=np.arange(0, XMAX, step)
         )
plt.xlabel("relative $\\beta$ beating")
plt.ylabel("Frequency")
plt.legend()
plt.savefig(f"recon_{TYPE_SORTING.lower()}.pdf")

plt.close()

# --------------------------------------------------------------------------------------------------
fig = plt.gcf()

plt.xlim(0,6)
plt.ylim(0,4)
fig.set_size_inches(6, 4)
plt.scatter(summ["BBEATX_AFTER"]/summ["BBEATX"], summ["CORRX_AFTER"]/summ["CORRX"],
            s=2)
plt.xlabel("BBEATX improvement")
plt.ylabel("CORRX improvement")

plt.savefig(f"correlation_{TYPE_SORTING.lower()}.pdf")

# print("=== Dif ================================")
# 
# summ = tfs.read("summ_diff.tfs")
# 
# data = summ["CORR_AFTER"]/summ["CORR"]
# 
# plt.close()
# 
# data_over = [x for x in data if x >= 1]
# print("CORR > 0")
# print(len(data_over)/len(data))
# 
# data_under = [x for x in data if x < 1]
# print("CORR < 0")
# print(len(data_under)/len(data)) 
# 
# plt.hist(data_under,
#          label="CORR",
#          color="#5050F0",
#          edgecolor="#000080",
#          histtype="stepfilled",
#          bins=np.arange(0, 1.1, STEP))
# 
# plt.hist(data_over,
#          label="CORR",
#          color="#F05050",
#          edgecolor="#800000",
#          histtype="stepfilled",
#          bins=np.arange(1, 8, STEP))
# plt.legend()
# plt.savefig("diff.pdf")
# 
# print(summ["CORR_AFTER"]/summ["CORR"])
# print(summ["BBEAT_AFTER"]/summ["BBEAT"])
