import psutil
import numpy as np
from subprocess import call
import time

parameters = [0.5, 1.0, 1.5]
while len(parameters) > 0:
  cores_usage = psutil.cpu_percent(interval=2, percpu=True)
  free_cores = np.where(np.asarray(cores_usage) < 50)
  free_cores = len(free_cores[0])
  if free_cores > 10:
      rc = call("./run_script.sh %s" % str(parameters[0]), shell = True)
      parameters.pop(0)
      time.sleep(10)

