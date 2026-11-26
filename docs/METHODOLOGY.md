# Historical methods and limitations

The provisional 67-cell CWatM catchment drains to global 0.5-degree row 110, column 505 (34.75 N, 72.75 E). The modeled basin area is approximately 171,171 km², versus 162,393 km² in the WAPDA station metadata. The RGI 6.2 screening intersects regions 13, 14 and 15 with this catchment and produced 13,952 candidates. The catchment and glacier release cells are not independently validated physical watershed boundaries.

ERA5 hourly instant variables are averaged on UTC days (temperature extrema and wind magnitude are evaluated before daily averaging); accumulated variables use hours 01 through the following 00. Specific humidity is calculated from pressure and dewpoint. Native 0.25-degree forcing is area-overlap remapped onto the 14 by 21 regional CWatM grid. Tiny negative shortwave radiation roundoff is screened using the historical 10 J/m² tolerance. The source dataset, geographic extent, period, variables and conversion choices are encoded in the climate modules.

Glacier production uses ERA5-calibrated OGGM 1.6.3 level-3 glacier directories and fixed model-flowline geometry. Daily climate variability is reconciled to the OGGM monthly climate and the monthly calibrated melt and phase totals; melt is partitioned into snow, firn and ice. This produces potential liquid-water supply, not independently measured daily glacier discharge. Snow/firn start empty in 2005; the historical method does not explicitly represent dynamic geometry, refreezing or delayed glacier water release. Seven failed reference glaciers are documented in the configuration and omitted (approximately 2.322 km² in-basin overlap; 99.9871% area coverage). The export scales by in-basin polygon overlap, capped at unity, assigns discharge to the largest-overlap candidate cell, and distributes cover according to polygon intersections. A substantial part of exported water is associated with provisional routing/membership flags. Daily melt and rain are regional 14×21 arrays, while cover is global; this avoids the proven CWatM regional reader mis-slice. The unit for daily glacier water inputs is m³/day.

CWatM uses a continuous 2005 spin-up, calibration 2006–2012, and evaluation 2013–2019. The selected precipitation multiplier is 0.70 (precipitation_coversion 60.48), IceMeltCoef 0.0035 and manningsN 2.0. Other settings remain unchanged. The precipitation adjustment is empirical and not a separately validated ERA5 correction; OGGM retains its ERA5 calibration. CWatM IceMelt represents additional summer snowmelt, not independently modeled glacier ice melt. The calibration search reproduces the historical tested cases and prioritizes annual peak magnitude error, not a claim of global optimality.

## Reference performance

| Period | NSE | KGE | Bias (%) | Mean absolute relative annual peak error (%) |
|---|---:|---:|---:|---:|
| Calibration 2006–12 | 0.828220 | 0.880978 | -1.110975 | 19.774970 |
| Evaluation 2013–19 | 0.818682 | 0.892321 | -0.723442 | 10.489692 |
| Full 2006–19 | 0.823338 | 0.887267 | -0.916742 | 15.132331 |

Annual peak timing error averages 11.93 days. The model produces 153 days above 10,000 m³/s versus 69 observed in the full period; these are exceedance days, not independent events. May–June low bias and March–April high bias remain. The evaluation years were not used to choose the final parameters, but their uncalibrated performance had previously been inspected. The historical fixed-geometry method must not simply be extrapolated to 2100. Future projections require consistent climate forcing and evolving glacier geometry, with uncertainty in absolute flood frequencies reported explicitly.
