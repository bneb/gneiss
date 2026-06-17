import json

# compute MW combination from the first few epochs of WTZR
c = 299792458.0
f1 = 1575.42e6
f2 = 1227.60e6
lam1 = c / f1
lam2 = c / f2

osb_L1 = 0.88372 * 1e-9 * c # G08 L1C/L1W
osb_L2 = 1.47093 * 1e-9 * c # G08 L2C/L2W
osb_P1 = 0.6038 * 1e-9 * c  # G08 C1C
osb_P2 = 1.2970 * 1e-9 * c  # G08 C2C
osb_P1W = 0.0
osb_P2W = 0.0

def mw(L1, L2, P1, P2):
    L1 = L1 * lam1
    L2 = L2 * lam2
    return (f1 * L1 - f2 * L2) / (f1 - f2) / (c / (f1 - f2)) - (f1 * P1 + f2 * P2) / (f1 + f2) / (c / (f1 - f2))

# Wait, I don't have python reading RNX easily. 
# Let me write a Rust script inside gneiss to compute MW!
