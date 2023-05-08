
#
# import cProfile, pstats, io
# # from pstats import SortKey
# pr = cProfile.Profile()
# pr.enable()
# # ... do something ...
#
# import compare_projectors as mod
# stats.sort_stats('cumulative')
# mod.run()
#
# pr.disable()
# s = io.StringIO()
# sortby = SortKey.CUMULATIVE
# ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
# ps.print_stats()
# print(s.getvalue())
#
#

### do this: from https://docs.python.org/3/library/profile.html#module-cProfile and https://pymotw.com/2/profile/

# in a terminal:

python -m cProfile -o projstats compare_projectors.py

# Read all 5 stats files into a single object
stats = pstats.Stats('projstats')

# Sort the statistics by the cumulative time spent in the function
stats.sort_stats('cumulative')

stats.print_stats(10)