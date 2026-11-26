"""Fixed-geometry OGGM monthly-constrained daily melt and rain disaggregation.

The monthly OGGM quantities are imposed, not independently validated. Daily
ERA5 anomalies and positive-degree-day shares determine intramonth timing.
The snow and five firn cohorts partition the melt; no refreezing or routing.
"""
from __future__ import annotations
import calendar, math
from datetime import date, timedelta
import numpy as np
from .common import START, END, DAYS as DATE_LIST
DAYS=len(DATE_LIST)
MONTHLY_ABS_TOL_M3=1e-6
DAILY_FIELDS=['date','snowfall_m3','rain_m3','snowmelt_m3','firnmelt_m3','icemelt_m3','total_melt_m3','surface_liquid_m3','snow_firn_storage_m3_end','firn_to_ice_m3']
class RecoveryError(RuntimeError): pass
class DiagnosticError(RuntimeError): pass

def glacier_bin_geometry(flowlines: list) -> tuple[np.ndarray, np.ndarray]:
    """Return heights and actual ice-covered areas of model flowline bins.

    A flowline bin with zero thickness is not part of the glacierized area.
    bin_area_m2 is an OGGM flowline property based on physical width/spacing.
    """
    hs, areas = [], []
    for fl in flowlines:
        h = np.asarray(fl.surface_h, dtype=np.float64)
        thick = np.asarray(fl.thick, dtype=np.float64)
        a = np.asarray(fl.bin_area_m2, dtype=np.float64)
        if h.ndim != 1 or h.shape != thick.shape or h.shape != a.shape:
            raise ValueError('Inconsistent OGGM flowline bin geometry')
        if not np.all(np.isfinite(h)) or not np.all(np.isfinite(thick)) or not np.all(np.isfinite(a)):
            raise ValueError('Nonfinite OGGM flowline geometry')
        if np.any(a < 0) or np.any(thick < 0):
            raise ValueError('Negative bin area or ice thickness')
        use = (thick > 0) & (a > 0)
        hs.extend(h[use])
        areas.extend(a[use])
    heights = np.asarray(hs, dtype=np.float64)
    area_m2 = np.asarray(areas, dtype=np.float64)
    if heights.size == 0 or area_m2.sum() <= 0:
        raise ValueError('No ice-covered flowline bins found')
    return heights, area_m2

def partition_precip(prcp, temp, t_solid, t_liq):
    """OGGM's linear temperature-dependent snowfall fraction; units kg m-2 day-1."""
    if not (np.isfinite(t_solid) and np.isfinite(t_liq) and t_liq > t_solid):
        raise DiagnosticError('Invalid solid/liquid thresholds')
    frac = np.clip((t_liq - temp) / (t_liq - t_solid), 0, 1)
    return prcp * frac, prcp * (1 - frac)

def climate_month(model, year: int, month: int, month_rows, heights: np.ndarray,
                  date_to_floatyear):
    """Construct daily fields matching the calibrated OGGM monthly background.

    This function does not itself force melt
    to agree with OGGM. Units: T degC; snow/rain kg m-2 day-1; PDD K day.
    """
    n = calendar.monthrange(year, month)[1]
    if len(month_rows) != n or any(d.year != year or d.month != month for d, _, _ in month_rows):
        raise DiagnosticError(f'Incomplete or mixed month {year}-{month:02d}')
    indexes = np.flatnonzero((np.asarray(model.years) == year) & (np.asarray(model.months) == month))
    if len(indexes) != 1:
        raise DiagnosticError(f'OGGM lacks unique monthly record {year}-{month:02d}')
    ix = int(indexes[0])
    t_raw = np.array([x[1] for x in month_rows], dtype=float)
    p_raw = np.array([x[2] for x in month_rows], dtype=float)
    p_month_target = float(model.prcp[ix])
    if not np.isfinite(p_month_target) or p_month_target < 0:
        raise DiagnosticError('Invalid OGGM calibrated monthly precipitation')
    p_sum = float(p_raw.sum())
    if p_sum == 0 and p_month_target > 1e-9:
        raise DiagnosticError(f'Cannot scale dry daily month to wet OGGM monthly {year}-{month:02d}')
    p_factor = p_month_target / p_sum if p_sum > 0 else 0.0
    t_month = float(model.temp[ix])
    lapse = float(model.grad[ix])
    if not np.isfinite([t_month, lapse, model.ref_hgt, model.t_melt]).all():
        raise DiagnosticError('Invalid OGGM calibrated monthly temperature')
    bins_t = (t_month + t_raw - t_raw.mean())[:, None] + (heights - model.ref_hgt)[None, :] * lapse
    bins_p = (p_raw * p_factor)[:, None] * np.ones((1, heights.size))
    snow, rain = partition_precip(bins_p, bins_t, float(model.t_solid), float(model.t_liq))
    pdd = np.maximum(bins_t - float(model.t_melt), 0.0)
    og_t, og_tmelt, og_p, og_snow = model.get_monthly_climate(heights, year=date_to_floatyear(year, month))
    if not np.allclose(bins_t.mean(axis=0), np.asarray(og_t), atol=1e-8, rtol=0):
        raise DiagnosticError('Daily temperature mean not consistent with OGGM')
    if not np.allclose(bins_p.sum(axis=0), np.asarray(og_p), atol=1e-6, rtol=1e-11):
        raise DiagnosticError('Daily precipitation does not sum to OGGM monthly')
    if not np.allclose(snow + rain, bins_p, atol=1e-10, rtol=0):
        raise DiagnosticError('Daily precipitation phase balance failure')
    if not np.isfinite(np.asarray(og_tmelt, dtype=float)).all():
        raise DiagnosticError('Nonfinite OGGM monthly melt')
    return {'days': [x[0] for x in month_rows], 'snow': snow, 'rain': rain,
            'pdd': pdd, 'monthly_melt_ref_per_bin': np.asarray(og_tmelt, dtype=float) * float(model.monthly_melt_f),
            'monthly_rain_ref_per_bin': np.asarray(og_p, dtype=float) - np.asarray(og_snow, dtype=float),
            'monthly_snow_ref_per_bin': np.asarray(og_snow, dtype=float),
            'monthly_precip_ref_per_bin': np.asarray(og_p, dtype=float),
            'precip_scaling': p_factor}

def checked_volume(values, areas):
    """m3 water equivalent from kg/m2 and m2; stable bin accumulation."""
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1 or arr.shape != areas.shape or not np.all(np.isfinite(arr)):
        raise RecoveryError('Invalid water-equivalent bin array')
    return math.fsum(float(x) * (float(a) / 1000.0) for x, a in zip(arr, areas))

def validate_fields(field, areas):
    days = field['days']
    count, bins = len(days), len(areas)
    if count < 28 or count > 31:
        raise RecoveryError('Unexpected month length')
    for key in ('pdd', 'snow', 'rain'):
        x = np.asarray(field[key], dtype=float)
        if x.shape != (count, bins) or not np.isfinite(x).all() or x.min() < -1e-9:
            raise RecoveryError(f'Invalid {key} matrix')
    for key in ('monthly_melt_ref_per_bin', 'monthly_rain_ref_per_bin',
                'monthly_snow_ref_per_bin', 'monthly_precip_ref_per_bin'):
        x = np.asarray(field[key], dtype=float)
        if x.shape != (bins,) or not np.isfinite(x).all() or x.min() < -1e-8:
            raise RecoveryError(f'Invalid {key} vector')
    return count, bins

def split_daily_month(field, areas):
    """Constrain month melt & rain at EACH elevation bin, retaining daily PDD/precip timing.

    Month precipitation is rescaled by band only to repair numerical differences
    from the calibrated OGGM precipitation total. An additive change in the
    day-wise rain fraction (clipped to [0,1]) reconciles the different daily
    and monthly temperature-dependent phase assumptions. All such changes are
    flagged: forcing monthly phase does not physically validate daily phase.
    """
    n, b = validate_fields(field, areas)
    pdd = np.maximum(np.asarray(field['pdd'], dtype=float), 0.0)
    unscaled_snow = np.maximum(np.asarray(field['snow'], dtype=float), 0.0)
    unscaled_rain = np.maximum(np.asarray(field['rain'], dtype=float), 0.0)
    unscaled_precip = unscaled_snow + unscaled_rain
    melt_ref = np.maximum(np.asarray(field['monthly_melt_ref_per_bin'], dtype=float), 0.0)
    rain_ref = np.maximum(np.asarray(field['monthly_rain_ref_per_bin'], dtype=float), 0.0)
    snow_ref = np.maximum(np.asarray(field['monthly_snow_ref_per_bin'], dtype=float), 0.0)
    precip_ref = np.maximum(np.asarray(field['monthly_precip_ref_per_bin'], dtype=float), 0.0)
    if abs(checked_volume(snow_ref + rain_ref - precip_ref, areas)) > MONTHLY_ABS_TOL_M3:
        raise RecoveryError('OGGM monthly snow/rain reference does not conserve precipitation')
    # Also check per-bin reference consistency; positive and negative mismatch
    # between bins must not be allowed to cancel in the volume sum.
    if checked_volume(np.abs(snow_ref + rain_ref - precip_ref), areas) > MONTHLY_ABS_TOL_M3:
        raise RecoveryError('OGGM monthly phase reference inconsistent by elevation bin')
    if checked_volume(np.maximum(rain_ref - precip_ref, 0), areas) > MONTHLY_ABS_TOL_M3:
        raise RecoveryError('OGGM monthly rain exceeds total precipitation')

    melt = np.zeros_like(pdd)
    rain = np.zeros_like(pdd)
    snowfall = np.zeros_like(pdd)
    phase_offsets = []
    forced_zero_phase = forced_full_phase = 0
    phase_corrected_bins = 0
    for j in range(b):
        # PDD support must exist whenever OGGM specifies positive melt.
        mref = float(melt_ref[j]); pddsum = math.fsum(map(float, pdd[:, j]))
        if mref > 0:
            if pddsum <= 0:
                raise RecoveryError(f'Unallocatable OGGM monthly melt: no positive daily PDD in bin {j} on {field["days"][0]}')
            melt[:, j] = pdd[:, j] * (mref / pddsum)
            k = int(np.argmax(melt[:, j]))
            melt[k, j] += mref - math.fsum(map(float, melt[:, j]))
            if (melt[:, j] < -1e-8).any():
                raise RecoveryError('Negative allocated melt')

        ptarget = float(precip_ref[j]); rtarget = float(rain_ref[j])
        raw = unscaled_precip[:, j]
        rawtotal = math.fsum(map(float, raw))
        if ptarget > 0 and rawtotal <= 0:
            raise RecoveryError(f'Unallocatable OGGM monthly precipitation: no wet daily forcing in bin {j}')
        if ptarget == 0:
            if rtarget > 1e-10:
                raise RecoveryError('Positive reference rain with zero reference precipitation')
            continue
        if rtarget > ptarget:
            if checked_volume(np.array([max(rtarget-ptarget,0.)]), np.array([areas[j]])) > MONTHLY_ABS_TOL_M3:
                raise RecoveryError('Monthly rain greater than precipitation in bin')
            rtarget = ptarget
        daily_total = raw * (ptarget / rawtotal)
        k = int(np.argmax(daily_total))
        daily_total[k] += ptarget - math.fsum(map(float, daily_total))
        frac = np.divide(unscaled_rain[:, j], raw, out=np.zeros(n), where=raw > 0)
        frac = np.clip(frac, 0.0, 1.0)
        raw_scaled_rain = math.fsum(float(x) for x in unscaled_rain[:, j]) * (ptarget / rawtotal)
        if abs(raw_scaled_rain - rtarget) * (float(areas[j]) / 1000.0) > MONTHLY_ABS_TOL_M3:
            phase_corrected_bins += 1
        if rtarget <= 0:
            rain[:, j] = 0.0
            phase_offsets.append(-float(np.max(frac)))
            if raw_scaled_rain > 1e-12: forced_zero_phase += 1
        elif rtarget >= ptarget:
            rain[:, j] = daily_total
            wet = daily_total > 0
            phase_offsets.append(float(1.0 - np.min(frac[wet])))
            if ptarget - raw_scaled_rain > 1e-12: forced_full_phase += 1
        else:
            # Monotone constrained phase projection: exact monthly rain while
            # keeping every day's precipitation total and daily ordering fixed.
            lo, hi = -1.0, 1.0
            for _ in range(75):
                mid = 0.5 * (lo + hi)
                trial = float(np.dot(daily_total, np.clip(frac + mid, 0.0, 1.0)))
                if trial < rtarget: lo = mid
                else: hi = mid
            offset = 0.5 * (lo + hi)
            phase_offsets.append(offset)
            rain[:, j] = daily_total * np.clip(frac + offset, 0.0, 1.0)
            # Tiny double-precision residual must be corrected without
            # changing any daily precipitation or creating negative rain.
            correction = rtarget - math.fsum(map(float, rain[:, j]))
            if correction >= 0:
                k = int(np.argmax(daily_total - rain[:, j]))
                if correction > daily_total[k] - rain[k, j] + 1e-8:
                    raise RecoveryError('Unable to close rain phase deficit')
            else:
                k = int(np.argmax(rain[:, j]))
                if -correction > rain[k, j] + 1e-8:
                    raise RecoveryError('Unable to close rain phase surplus')
            rain[k, j] += correction
        snowfall[:, j] = daily_total - rain[:, j]
        if min(float(np.min(snowfall[:, j])), float(np.min(rain[:, j]))) < -1e-8:
            raise RecoveryError('Phase projection generated negative rainfall/snowfall')

    melt_m3 = checked_volume(melt_ref, areas)
    rain_m3 = checked_volume(rain_ref, areas)
    snow_m3 = checked_volume(snow_ref, areas)
    precip_m3 = checked_volume(precip_ref, areas)
    original_daily_rain_m3 = math.fsum(checked_volume(unscaled_rain[i], areas) for i in range(n))
    # Close at the *glacier-wide monthly* level, as requested, in absolute m3.
    for label, field_daily, wanted in (('melt', melt, melt_m3),
                                       ('rain', rain, rain_m3),
                                       ('snow', snowfall, snow_m3),
                                       ('precip', rain + snowfall, precip_m3)):
        obtained = math.fsum(checked_volume(field_daily[i], areas) for i in range(n))
        error = obtained - wanted
        if abs(error) > MONTHLY_ABS_TOL_M3:
            raise RecoveryError(f'{label} absolute monthly closure {error:.12g} m3 exceeds 1e-6 m3 on {field["days"][0]}')
    return melt, rain, snowfall, {
        'month': field['days'][0].strftime('%Y-%m'),
        'reference_melt_m3': melt_m3,
        'reference_rain_m3': rain_m3,
        'reference_snow_m3': snow_m3,
        'reference_precip_m3': precip_m3,
        'original_daily_rain_m3': original_daily_rain_m3,
        'phase_reconciliation_rain_difference_m3': rain_m3 - original_daily_rain_m3,
        'maximum_absolute_phase_fraction_offset': max(map(abs, phase_offsets), default=0.0),
        'rain_zero_reference_bins': forced_zero_phase,
        'rain_full_reference_bins': forced_full_phase,
        'rain_phase_corrected_bins': phase_corrected_bins,
    }

def glacier_daily(fields, areas):
    """Compute constrained potential and snow/firn/ice bookkeeping, 2005--2019."""
    areas = np.asarray(areas, dtype=float)
    if areas.ndim != 1 or not areas.size or not np.isfinite(areas).all() or (areas <= 0).any():
        raise RecoveryError('Invalid fixed OGGM elevation-band glacier geometry')
    if len(fields) != 180:
        raise RecoveryError('Expected 180 months of 2005--2019 OGGM input')
    snowpack = np.zeros(len(areas), dtype=float)
    firn = np.zeros((5, len(areas)), dtype=float)
    daily, monthly = [], []
    prev = START
    for f in fields:
        melt, rain, snowfall, diag = split_daily_month(f, areas)
        month_daily = []
        for i, d in enumerate(f['days']):
            if d != prev:
                raise RecoveryError(f'Daily date gap/duplicate; expected {prev}, got {d}')
            prev += timedelta(days=1)
            m = melt[i]
            s = snowfall[i]
            r = rain[i]
            before = checked_volume(snowpack + firn.sum(axis=0), areas)
            available_snow = snowpack + s
            snowmelt = np.minimum(available_snow, m)
            snowpack = available_snow - snowmelt
            remainder = m - snowmelt
            firnmelt = np.zeros(len(areas), dtype=float)
            for age in range(5):
                take = np.minimum(firn[age], remainder)
                firn[age] -= take
                firnmelt += take
                remainder -= take
            icemelt = np.maximum(remainder, 0.0)
            snowvol, rainvol = checked_volume(s, areas), checked_volume(r, areas)
            sv, fv, iv = (checked_volume(x, areas) for x in (snowmelt, firnmelt, icemelt))
            to_ice = 0.0
            if d.month == 12 and d.day == 31:
                to_ice = checked_volume(firn[-1], areas)
                firn[1:] = firn[:-1].copy()
                firn[0] = snowpack
                snowpack = np.zeros_like(snowpack)
            after = checked_volume(snowpack + firn.sum(axis=0), areas)
            water_err = (before + snowvol) - (sv + fv + to_ice + after)
            # Exact cancellation of huge stocks can exceed 1e-6 m3; keep a
            # separate scale-aware numerical *reservoir* condition. Monthly
            # melt/rain/precip budgets use the strict 1e-6 m3 absolute gate.
            if abs(water_err) > max(1e-6, 1e-11 * max(1., abs(before), abs(after), snowvol)):
                raise RecoveryError(f'Snow/firn accounting closure failed on {d}: {water_err} m3')
            total_melt = sv + fv + iv
            row = {'date': d.isoformat(), 'snowfall_m3': snowvol, 'rain_m3': rainvol,
                   'snowmelt_m3': sv, 'firnmelt_m3': fv, 'icemelt_m3': iv,
                   'total_melt_m3': total_melt,
                   'surface_liquid_m3': total_melt + rainvol,
                   'snow_firn_storage_m3_end': after, 'firn_to_ice_m3': to_ice}
            if any(not math.isfinite(v) or v < -1e-8 for k, v in row.items() if k != 'date'):
                raise RecoveryError(f'Negative/nonfinite daily flux/storage on {d}')
            month_daily.append(row)
            daily.append(row)
        residuals = {}
        for metric, ref in (('total_melt_m3', diag['reference_melt_m3']),
                            ('rain_m3', diag['reference_rain_m3']),
                            ('snowfall_m3', diag['reference_snow_m3'])):
            actual = math.fsum(float(row[metric]) for row in month_daily)
            residuals[metric] = actual - ref
            if abs(residuals[metric]) > MONTHLY_ABS_TOL_M3:
                raise RecoveryError(f'Monthly {metric} output closure {residuals[metric]:.12g} m3 exceeds 1e-6 m3 on {diag["month"]}')
        diag['monthly_output_residual_m3'] = residuals
        monthly.append(diag)
    if prev != END + timedelta(days=1) or len(daily) != DAYS:
        raise RecoveryError('Daily output does not cover every day 2005--2019')
    metrics = ('total_melt_m3', 'rain_m3', 'snowfall_m3')
    residual_max = {metric: max(abs(x['monthly_output_residual_m3'][metric]) for x in monthly) for metric in metrics}
    return daily, monthly, {
        'monthly_max_abs_residual_m3': residual_max,
        'monthly_absolute_tolerance_m3': MONTHLY_ABS_TOL_M3,
        'total_reference_melt_m3': math.fsum(x['reference_melt_m3'] for x in monthly),
        'total_reference_rain_m3': math.fsum(x['reference_rain_m3'] for x in monthly),
        'total_disaggregated_melt_m3': math.fsum(x['total_melt_m3'] for x in daily),
        'total_disaggregated_rain_m3': math.fsum(x['rain_m3'] for x in daily),
        'maximum_phase_fraction_offset': max(x['maximum_absolute_phase_fraction_offset'] for x in monthly),
        'monthly_rain_phase_adjustment_m3_total_abs': math.fsum(abs(x['phase_reconciliation_rain_difference_m3']) for x in monthly),
        'daily_days': DAYS,
    }
