# Example calculation using pwr package
library(pwr)

# Power to detect your observed MMP2 pattern (4/4 vs 0/5)
pwr.2p.test(h = ES.h(1.0, 0.0),
            n = 4, # Use smaller group size
            sig.level = 0.05)

# Power for moderate effects like CHI3L1 (3/4 vs 1/5)
pwr.2p.test(h = ES.h(0.75, 0.20),
            n = 4,
            sig.level = 0.05)

# Calculate required sample size for 80% power to detect moderate effects
pwr.2p.test(h = ES.h(0.75, 0.20),
            sig.level = 0.05,
            power = 0.8)
# This will show you need ~12-15 per group for adequate power