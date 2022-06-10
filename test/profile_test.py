import os
import sys
import pstats
from pstats import SortKey

to_profile = sys.argv[1]
if len(sys.argv) > 2:
    percentage = int(sys.argv[2])
else:
    percentage = 5

print(f"Profiling {to_profile}")
os.system(f"python -m cProfile -o profile.out {to_profile}")
os.system(f"gprof2dot -f pstats -n {percentage} profile.out | dot -Tpng -o profile.png && eog profile.png")

