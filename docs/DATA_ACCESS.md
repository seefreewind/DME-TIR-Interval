# Data access and source attribution

Obtain the original packages directly from the DRCR public-data portal, follow its registration and usage requirements, and read each package README and glossary before use:

- Protocol T: https://public.jaeb.org/drcrnet/stdy/download/206 → `data/raw/T/Protocol_T.zip`
- Protocol AC: https://public.jaeb.org/drcrnet/stdy/download/514 → `data/raw/AC/Protocol_AC.zip`

Verified archive SHA256 values:

| Trial | SHA256 |
| --- | --- |
| T | be6be58ec51e9d257ee6d6881fb4e22d605ff64d79483435ee003b0fefffe9c5 |
| AC | 49f4dac3b2a3324c5ebab7f7d287840b4102ed23aae110578b0769badf955a1a |

The code fails closed if a package hash or frozen derived-input hash differs. A revised official package requires an explicit adaptation audit; do not bypass the check. All packages, extracted files, participant identifiers, eye-level scores, model weights and individual predictions remain local and are excluded from Git and Zenodo. The MIT license applies to original project code, not the source datasets or third-party materials.

## Protocol T attribution required by its README

The source of the data is the DRCR Retina Network (Unique Federal Award Identification Number (FAIN) UG1EY014231). Protocol T: A Comparative Effectiveness Study of Intravitreal Aflibercept, Bevacizumab and Ranibizumab for Diabetic Macular Edema [Data file]. Retrieved from https://public.jaeb.org/drcrnet. The analyses content and conclusions presented herein are solely the responsibility of the authors and have not been reviewed or approved by the DRCR Retina Network.

## Protocol AC attribution required by its README

The source of the data is the DRCR Retina Network (Unique Federal Award Identification Number (FAIN) UG1EY014231). Protocol AC: Randomized Trial of Intravitreous Aflibercept versus Intravitreous Bevacizumab + Deferred Aflibercept for Treatment of Central-Involved Diabetic Macular Edema [Data file]. Retrieved from https://public.jaeb.org/drcrnet. The analyses content and conclusions presented herein are solely the responsibility of the authors and have not been reviewed or approved by the DRCR Retina Network or Jaeb Center for Health Research.

The dataset-source award number above is attribution of the original dataset and does not indicate funding of this secondary analysis.
