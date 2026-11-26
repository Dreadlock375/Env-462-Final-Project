"""Fast numerical regression checks using archived reference tables."""
import datetime as dt
import re
import tempfile
from unittest.mock import patch
from pathlib import Path
import unittest
import numpy as np
from karakoram import basin, coupling, glaciers, hydrology, observations, prepare_climate, setup_inputs
from karakoram.common import ROOT, CONFIG, read_csv, read_json

class ReferenceTests(unittest.TestCase):
    @unittest.skipUnless(observations.SOURCE.is_file(), "Original WAPDA data required for observed-metric regression")
    def test_frozen_metrics(self):
        rows=read_csv(ROOT/'reference/validation_daily.csv')
        dates=[dt.date.fromisoformat(r['date']) for r in rows]
        records,_=observations.read_daily_csv(observations.SOURCE)
        s=np.array([float(r['calibrated_m3s']) for r in rows]); o=np.array([records[d]['discharge_m3s'] for d in dates])
        hydrology.verify_reference_observations(o)
        expected=read_json(ROOT/'reference/validation_report.json')['periods']
        n=2557
        for key,sl in [('calibration_2006_2012',slice(0,n)),('validation_2013_2019',slice(n,None)),('full_2006_2019',slice(None))]:
            actual=hydrology.metrics(s[sl],o[sl],dates[sl])
            for metric in ('nse','kge','bias_pct','annual_peak_abs_relative_error_pct','annual_peak_mean_abs_timing_error_days'):
                self.assertAlmostEqual(actual[metric],expected[key]['calibrated'][metric],places=9)
    def test_spatial_reference(self):
        jobs=read_csv(ROOT/'reference/glacier_worklist.csv'); self.assertEqual(len(jobs),13952)
        self.assertEqual(len(set(j['rgi_id'] for j in jobs)),13952)
        cells=read_csv(ROOT/'reference/basin_cells.csv');self.assertEqual(len(cells),67)
        overlap=read_csv(ROOT/'reference/glacier_overlap.csv');self.assertEqual(len(overlap),14397)
        self.assertAlmostEqual(sum(float(j['polygon_overlap_besham_km2']) for j in jobs),18056.292811,places=5)
    def test_original_worklist(self):
        candidates=read_csv(ROOT/'reference/glacier_candidates.csv')
        overlap=read_csv(ROOT/'reference/glacier_overlap.csv')
        cells=read_csv(ROOT/'reference/basin_cells.csv')
        mask=np.zeros((360,720),dtype=bool)
        for r in cells:mask[int(r['row']),int(r['col'])]=True
        screen,counts=basin.verify_candidate_tables(candidates,overlap,mask)
        jobs=basin.make_jobs(candidates,screen)
        self.assertEqual([{k:str(v) for k,v in j.items()} for j in jobs],read_csv(ROOT/'reference/glacier_worklist.csv'))
    @unittest.skipUnless(observations.SOURCE.is_file(), "Original WAPDA data required for observation regression")
    def test_raw_observations(self):
        records,_=observations.read_daily_csv(observations.SOURCE)
        rows=read_csv(ROOT/'reference/validation_daily.csv')
        self.assertEqual(len(records),7305)
        values=[records[dt.date.fromisoformat(r['date'])]['discharge_m3s'] for r in rows]
        hydrology.verify_reference_observations(values)
        values[0]+=1.0
        with self.assertRaisesRegex(RuntimeError,'Reference observations differ'):
            hydrology.verify_reference_observations(values)
    def test_routing_graph(self):
        # Synthetic grid: all three upstream cells reach the bottom-right sink.
        ldd=np.array([[3.,2.],[6.,5.]])
        area=np.ones((2,2));lat=np.array([1.,0.]);lon=np.array([0.,1.])
        net=basin.build_network(ldd,area,lat,lon)
        self.assertTrue(basin.upstream_mask(net.downstream,net.valid,3).all())
        self.assertEqual(net.upstream_area_m2[3],4.)
    def test_daily_climate(self):prepare_climate.self_test()
    def test_coupling(self):coupling.self_test()
    def test_glacier(self):glaciers.self_test()
    def test_settings(self):
        text=hydrology.settings(ROOT/'results/test',hydrology.LAST,hydrology.FROZEN)
        self.assertNotIn('@ROOT@',text);self.assertNotIn('@OUTPUT@',text)
        self.assertIn('precipitation_coversion = 60.48',text)
        self.assertRegex(text,r'IceMeltCoef\s*=\s*0\.0035')
        self.assertIn('manningsN = 2.0',text)
        # A portable template must not contain machine-specific absolute paths.
        # A rendered configuration, in contrast, must use this project's paths.
        template=(ROOT/'config/cwatm.ini.template').read_text()
        self.assertIn('@ROOT@',template)
        self.assertIn('@OUTPUT@',template)
        self.assertNotRegex(template,r'(?m)^\s*[^#;\[\n][^=\n]*=\s*/')
        absolute_paths=[]
        for line in text.splitlines():
            match=re.match(r'^\s*[^#;\[\n][^=\n]*=\s*(/.*)$',line)
            if match:
                path=Path(match.group(1).strip())
                self.assertTrue(path.is_relative_to(ROOT),f'Path outside project: {path}')
                absolute_paths.append(path)
        self.assertIn(ROOT/'results/test',absolute_paths)
        self.assertIn(ROOT/'data/processed/glacier_forcing/melt_on_glacier_2005_2019.nc',absolute_paths)

    def test_independent_input_copy(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory);root=base/'new';root.mkdir()
            source=base/'old';source.mkdir();(source/'input.txt').write_text('original data')
            with patch.object(setup_inputs,'ROOT',root):
                setup_inputs.copy_tree(source,root/'input')
                target=root/'input/input.txt'
                self.assertFalse(target.is_symlink())
                self.assertEqual(target.read_text(),'original data')
                setup_inputs.copy_file(source/'input.txt',target)
                target.write_text('different')
                with self.assertRaises(RuntimeError):setup_inputs.copy_file(source/'input.txt',target)
                (root/'link').symlink_to(source,target_is_directory=True)
                with self.assertRaises(RuntimeError):setup_inputs.copy_file(source/'input.txt',root/'link/input.txt')

    def test_local_oggm_acquisition_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);target=root/'data/raw/rgi'
            target.mkdir(parents=True)
            for region in CONFIG['rgi_regions']:
                for ext in ('.shp','.shx','.dbf','.prj'):
                    (target/(region+'_rgi62_test'+ext)).write_text('fixture')
            with patch.object(setup_inputs,'ROOT',root):
                setup_inputs.prepare_oggm('basin',allow_download=False)
                self.assertIn('has_internet = False',(root/'.runtime-home/.oggm_config').read_text())
                (target/'13_rgi62_test.shp').unlink()
                with self.assertRaisesRegex(RuntimeError,'explicitly allow acquisition'):
                    setup_inputs.prepare_oggm('basin',allow_download=False)
                self.assertIn('has_internet = False',(root/'.runtime-home/.oggm_config').read_text())

    def test_archive_groups_and_offline_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);worklist=root/'data/processed/basin/worklist.csv'
            worklist.parent.mkdir(parents=True)
            worklist.write_bytes((ROOT/'reference/glacier_worklist.csv').read_bytes())
            with patch.object(setup_inputs,'ROOT',root):
                archives=setup_inputs.glacier_archives(True)
                self.assertEqual(len(archives),35)
                for url,path in archives.items():
                    self.assertTrue(path.is_relative_to(root/'.runtime-home/OGGM/download_cache'))
                    self.assertTrue(url.endswith(path.name))
                    path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(b'fixture')
                with patch.object(setup_inputs,'auxiliary_present',return_value=True):
                    setup_inputs.prepare_oggm('glaciers',True,False)
                    next(iter(archives.values())).unlink()
                    with self.assertRaisesRegex(RuntimeError,'explicitly allow acquisition'):
                        setup_inputs.prepare_oggm('glaciers',True,False)
                self.assertIn('has_internet = False',(root/'.runtime-home/.oggm_config').read_text())

    def test_static_archive_verification(self):
        import zipfile
        from karakoram.common import sha
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);base=root/setup_inputs.STATIC_ROOT
            base.mkdir(parents=True)
            name=setup_inputs.STATIC_DIRS[0]
            archive=base/'test.zip';source=base/'fixture.nc';source.write_bytes(b'map')
            manifest={'static_sha256':{str(setup_inputs.STATIC_ROOT/name/'fixture.nc'):sha(source)}}
            with zipfile.ZipFile(archive,'w') as z:z.writestr(name+'/fixture.nc',b'map')
            with patch.object(setup_inputs,'ROOT',root):
                setup_inputs.extract_static_archive(archive,base/name,manifest)
                self.assertEqual(setup_inputs.static_subset(base/name,name,manifest),1)
                (base/name/'fixture.nc').write_bytes(b'changed')
                with self.assertRaisesRegex(RuntimeError,'Static map mismatch'):
                    setup_inputs.static_subset(base/name,name,manifest)
                with patch('urllib.request.urlopen',side_effect=AssertionError('No network expected')):
                    setup_inputs.verified_download('https://example.invalid/',source,sha(source))
                    with self.assertRaisesRegex(RuntimeError,'checksum differs'):
                        setup_inputs.verified_download('https://example.invalid/',source,'0'*64)
    def test_static_archive_rejects_traversal(self):
        import zipfile
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);archive=root/'test.zip'
            with zipfile.ZipFile(archive,'w') as z:z.writestr('../unsafe.nc',b'map')
            with patch.object(setup_inputs,'ROOT',root):
                with self.assertRaisesRegex(RuntimeError,'Unsafe ZIP member'):
                    setup_inputs.extract_static_archive(archive,root/'maps',{})
            self.assertFalse((root/'maps').exists())
    def test_captured_source_locks(self):
        manifest=read_json(ROOT/'config/source_manifest.json')
        self.assertEqual(manifest['cwatm_revision'],CONFIG['cwatm_revision'])
        self.assertEqual(len(manifest['static_sha256']),141)
        for name in ('oggm_env','cwatm_env'):
            lock=(ROOT/'config/locks'/f'{name}.conda.txt').read_text()
            self.assertIn('@EXPLICIT',lock)
            self.assertNotIn('/home/dreadlock/',lock)
            self.assertTrue(read_json(ROOT/'config/locks'/f'{name}.versions.json'))

if __name__=='__main__':unittest.main()
