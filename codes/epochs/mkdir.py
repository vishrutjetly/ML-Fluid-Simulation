import os

epochs = [1, 10, 100, 1000, 5000, 7000, 10000, 20000, 50000, 70000]

for e in epoch:
	os.mkdir("model_"+str(e))