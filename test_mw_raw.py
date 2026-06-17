import json

c = 299792458.0
f1 = 1575.42e6
f2 = 1227.60e6
lam1 = c / f1
lam2 = c / f2
lam_wl = c / (f1 - f2)

# From mw_test.rs G08
# Let's say we have raw values
# I don't have raw values printed. Let me write a tiny python script that parses rnx using simple regex just for the first epoch to get the raw numbers!

