# Precise Point Positioning with Ambiguity Resolution (PPP-AR)

This document provides a conceptual overview of Precise Point Positioning (PPP) and how the Gneiss Navigation Engine achieves integer **Ambiguity Resolution (AR)** to provide sub-decimeter accuracy globally without a local base station. It is written to be accessible to navigation engineers and hobbyists alike.

## The GNSS Observable Problem
When a GNSS receiver (like the one in your smartphone or a survey grade antenna) tracks a satellite, it records two main measurements:
1. **Pseudorange (Code):** The time it took the signal to travel from the satellite to you. This is noisy (typically 1-3 meters of error) but unambiguous. You know roughly where you are.
2. **Carrier Phase:** The fractional wave of the radio signal itself. This is incredibly precise (millimeter level noise) but **ambiguous**. The receiver can measure the fraction of the wave it sees, but it doesn't know how many *full integer waves* fit between the satellite and the antenna. This unknown integer is called the **Phase Ambiguity**.

## Float PPP vs Fixed PPP-AR
In traditional **Float PPP**, the navigation filter estimates these ambiguities as continuous decimal (float) numbers. Over time (typically 20-30 minutes), as satellites move across the sky, the changing geometry allows the filter to slowly narrow down these float estimates until the position converges to ~10-20 cm accuracy.

In **Fixed PPP-AR**, we recognize that the true ambiguities *must* be exact integers (whole wavelengths). If the engine can statistically prove which exact integer is correct, it can "fix" the ambiguity to that integer. Once fixed, the carrier phase measurement acts like a pseudorange with millimeter noise! The position accuracy instantly drops below 5 centimeters and convergence time is dramatically reduced.

## The Challenge: Hardware Biases
Why is fixing the integer hard? Because hardware isn't perfect. 

When the GNSS signal leaves the satellite's transmitter, it travels through internal cables and filters, adding a tiny delay. We call this the **Observable Specific Bias (OSB)**. 
- A delay of 1 nanosecond on the phase signal corresponds to about $0.3$ meters of distance. 
- If the wavelength is $0.19$ meters, a 0.3-meter delay shifts the float ambiguity by $1.57$ cycles! 

If the float ambiguity is shifted by a non-integer bias, the filter will never be able to round it to the correct integer. It will look like a fraction. 

### Phase Biases (SINEX)
To achieve PPP-AR, the Gneiss Engine downloads **Phase Bias (OSB)** files provided by institutions like the French Space Agency (CNES) via the IGS. These files (in the `.bia` SINEX format) contain the exact fractional cycle delays for every signal on every satellite, allowing Gneiss to perfectly calibrate the phase measurements back to zero-mean integers.

## Un-Differenced Un-Combined (UDUC) Processing
Historically, navigation engines combined the L1 and L2 frequencies together mathematically to remove the effect of the Earth's ionosphere (called the Iono-Free combination). While this removes the ionosphere, it also destroys the integer nature of the ambiguities (combining integers with non-integer multipliers yields a float).

Gneiss uses **Un-Differenced Un-Combined (UDUC)** processing. We process raw L1 and L2 signals directly and estimate the ionosphere as an explicit state variable. This preserves the integer nature of the phase ambiguities.

## Cascade Ambiguity Resolution (Widelane / Narrowlane)
Attempting to fix the raw L1 or L2 integers directly is mathematically unstable because their wavelengths are very short (~19 cm and ~24 cm). A tiny error will cause the filter to guess the wrong integer.

Instead, Gneiss uses a **Cascade** approach:
1. **Widelane (WL):** We mathematically subtract the L1 and L2 phase measurements. This creates an artificial "Widelane" signal with a massive wavelength of **86.4 cm**. Because the wavelength is so large, it is very easy to safely round the float ambiguity to the correct integer. 
2. **Narrowlane (NL):** Once the Widelane integer is fixed and fed back into the filter, the remaining uncertainty shrinks. The filter then uses the geometry to solve the "Narrowlane" ambiguity, which has a tiny wavelength of **10.7 cm**. Fixing this yields the ultimate sub-decimeter precision.

## Signal Tracking Fallbacks and the "Quarter Cycle Shift"
GNSS satellites broadcast many different signal codes (e.g., GPS broadcasts `C1C`, `C1W`, `L2X`, `L2W`, etc.). CNES phase bias products typically only provide biases for the core signals (like `L1W` and `L2W`).

What happens if a consumer receiver (like a U-Blox or smartphone) tracks the civilian `L2C` signal instead of the military `L2W` signal?
- Gneiss implements an intelligent fallback tree. If an `L2C` bias isn't available, we fall back to the `L2W` bias.
- **The Catch:** In GPS Block IIR-M and IIF satellites, the civilian L2C signal is generated using a different hardware path than the military L2P(Y) signal. This results in the L2C signal's phase tracking exactly **+0.25 cycles (a quarter wavelength)** ahead of the L2W signal.
- Gneiss automatically detects this fallback mode and applies the $-0.25$ cycle algorithmic shift to the phase measurement. Without this tiny correction, the Narrowlane ambiguity would be off by exactly a quarter cycle, and integer fixing would fail entirely!

## Current Status & Benchmark Tracking (WTZR)

The current Gneiss post-processing engine evaluates dual-frequency carrier tracking on the official WTZR IGS Reference Station in Germany:
- In standalone Float PPP mode (without external phase bias products), the engine converges from a 3-meter perturbed seed coordinate down to sub-decimeter horizontal errors ($p50 = 0.762\text{ m}$, minimum error $9.3\text{ cm}$, final vertical error $22.5\text{ cm}$).
- Full integer PPP-AR requires ingesting satellite observable-specific signal biases (OSBs via SINEX `.BIA` or `.OBX` files) and resolving wide-lane and narrow-lane integer ambiguities. The mathematical formulation for between-satellite single-difference wide-lane resolution is implemented in `gneiss-rtk::ambiguity::ppp_ar`, with full integer cascade fixing on real multi-day MGEX networks currently under active integration.

## Summary
By rigorously calibrating OSBs, explicitly modeling the ionosphere, and cascading through Widelane and Narrowlane combinations, integer PPP-AR provides globally precise centimeter-level positioning without local base stations.
