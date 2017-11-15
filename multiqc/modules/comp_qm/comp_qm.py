#usr/bin/env python

""" MultiQC module to parse output from Computational QM calculations e.g. Gaussian """

from __future__ import print_function
from collections import OrderedDict
import logging
import re
import json

from multiqc import config
from multiqc.plots import table, bargraph, linegraph, scatter, heatmap, beeswarm
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """Computational QM Analysis module"""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='computational_quantum', anchor='comp_qm',
        href="",
        info="Codes like Gaussian solve the Schrodinger equation, the aim here is to extract distance and charge measurements from the outputs")

        # Set up data structure
        self.comp_qm_data = {
            'mulliken': {},
            'bonds': {},
            'angles': {},
            'torsions': {}
        }

        # Find and load an  reports
        for f in self.find_log_files('comp_qm', filehandles=True):
            self.parse_comp_qm_log(f)

        # Filter to strip out ignored sample names
        self.comp_qm_data = self.ignore_samples(self.comp_qm_data)

        if len(self.comp_qm_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.comp_qm_data)))

        # Write parsed report data to a file
        self.write_data_file(self.comp_qm_data, 'multiqc_comp_qm')


        # Simple Table & Beeswarm for Mulliken charges
        if len(self.comp_qm_data['mulliken'])>0:
            self.comp_qm_mulliken_chart()

        # Simple Table and Beeswarm for Geometrical descriptors
        if len(self.comp_qm_data['bonds'])>0:
            self.comp_qm_geometry_descriptor_chart()
        #.. Commented out angles and dihedrals due to performance issues with highcharts - to work on
        # if len(self.comp_qm_data['angles'])>0:
            # self.comp_qm_geometry_descriptor_chart(geometry_descriptor="angles",min=0,max=180)
        # if len(self.comp_qm_data['torsions'])>0:
            # self.comp_qm_geometry_descriptor_chart(geometry_descriptor="torsions",min=-180.0,max=180.0)

    def parse_comp_qm_log(self, f):
        s_name = None

        s_name = self.clean_s_name(f['s_name'],f['root'])

        log.debug(str(s_name))
        log.debug("file contents")
        read_data=f['f'].read() # would be a bad idea for very large files
        log.debug(type(read_data))

        if s_name is not None:
            #. --- MATCH MULLIKEN BLOCKS---
            for match in re.finditer(r' Mulliken atomic charges.*?Sum of Mulliken atomic',read_data,re.DOTALL):
                pass # the last match is the one I want
            try:
                log.debug("Mulliken match information: %s", match)
                mulliken_charge_dictionary={} # dictionary to place parsed mulliken information, before placing in data source
                if match is not None:
                    log.debug(match.group(0)) # final match is what I want
                    all_charges=match.group(0).split('\n')[2:-1] # split into lines and remove unneccesary lines
                    for item in all_charges:
                        log.debug("Mulliken atomic number: %s", item.split()[0])
                        log.debug("Mulliken charge: %s",item.split()[2])
                        mulliken_charge_dictionary[str(item.split()[0])]=item.split()[2]

                    #. add data to the mulliken section
                    self.add_data_source(f, section='mulliken')
                    self.comp_qm_data['mulliken'][s_name] = mulliken_charge_dictionary
                    log.debug(str(self.comp_qm_data['mulliken'][s_name]))

            except:
                log.debug("Mulliken: No Data found.")

            #. --- MATCH DIST, ANGLE, DIHE ---
            for matchr in re.finditer(r' Optimized Parameters.*?GradGrad',read_data,re.DOTALL):
                pass # the last match is the one I want
            try:
                log.debug("Geometric Parameters match information: %s", match)
                #..  dictionaries to place parsed mulliken information, before placing in data source
                bond_dictionary={}
                angle_dictionary={}
                torsion_dictionary={}
                if matchr is not None:
                    log.debug(matchr.group(0)) # final match is what I want
                    all_vars=matchr.group(0).split('\n')[5:-2] # split into lines and remove unneccesary lines (-2 includes all)
                    for item in all_vars:
                        log.debug("Parameter name: %s",item.split()[1])
                        log.debug("Parameter value: %s", item.split()[3])
                        if "R" in item.split()[1]:
                            bond_dictionary[str(item.split()[1])]=item.split()[3]
                        elif "A" in item.split()[1]:
                            angle_dictionary[str(item.split()[1])]=item.split()[3]
                        elif "D" in item.split()[1]:
                            torsion_dictionary[str(item.split()[1])]=item.split()[3]
                        else:
                            log.debug("Failed the bond, angle, torsion check - unknown parameter")

                    #. add data to the bonds, angles, torsion section
                    self.add_data_source(f, section='bonds')
                    self.comp_qm_data['bonds'][s_name] = bond_dictionary
                    log.debug(str(self.comp_qm_data['bonds'][s_name]))
                    self.add_data_source(f, section='angles')
                    self.comp_qm_data['angles'][s_name] = angle_dictionary
                    log.debug(str(self.comp_qm_data['angles'][s_name]))
                    self.add_data_source(f, section='torsions')
                    self.comp_qm_data['torsions'][s_name] = torsion_dictionary
                    log.debug(str(self.comp_qm_data['torsions'][s_name]))
            except:
                log.debug("Parameters: No Data found.")
        s_name = None

    def comp_qm_mulliken_chart (self):
        """ Make the mulliken section table and plots """
        headers = OrderedDict()
        sample_keys=list(self.comp_qm_data['mulliken'].keys())

        for key in self.comp_qm_data['mulliken'][sample_keys[0]].keys(): # hack to get keys from first data set
            headers[key] = {
                'title': key,
                'description': 'Mulliken Charge',
                'suffix': '',
                'scale': 'Spectral', # 'RdBu',
                'dmin': -2.0,  # colour scale is set with dmin rather than min? see table.py  c_scale = mqc_colour.mqc_colour_scale(header['scale'], header['dmin'], header['dmax'])
                'dmax': 2.0,
                'ceiling': 2.0,
                'floor': -2.0,  # known issue with a negative range
                'format':'{:,.2f}',
                'shared_key':'mulliken_range'
            }

        log.debug(str(headers))


        # Config for the plot
        config = {
            'namespace': 'comp_qm',
            'id': 'comp_qm_mulliken',
            'title': 'Mulliken Charges',
        }
        self.add_section (
            name = 'Mulliken Charges',
            anchor = 'comp_qm_mulliken_table',
            plot = table.plot(self.comp_qm_data['mulliken'], headers, config)
        )
        #.. using dmin to ensure correct colour scales for negative ranges, note that min is still required for the beeswarm plot and that behaviour cancels out (?)
        #.. TLDR, the benefit of specifying a range is to keep all the graphs consistent. Maybe require an automated range check, rather than hardcoded. Try without if you like.
        #.. still seems strange, what is wrong. see table.py  c_scale = mqc_colour.mqc_colour_scale(header['scale'], header['dmin'], header['dmax'] and ../utils/mqc_colour.py
        for key in self.comp_qm_data['mulliken'][sample_keys[0]].keys(): # hack to get keys from first data set
            headers[key] = {
                'title': key,
                'description': 'Mulliken Charge',
                'suffix': '',
                'scale': 'Spectral', # 'RdBu',
                'min': -2.0,  # colour scale is set with dmin rather than min? see table.py  c_scale = mqc_colour.mqc_colour_scale(header['scale'], header['dmin'], header['dmax'])
                'max': 2.0,
                'ceiling': 2.0,
                'floor': -2.0,  # known issue with a negative range
                'format':'{:,.2f}',
                'shared_key':'mulliken_range'
            }

        self.add_section (
            name = 'Mulliken Charges: Beeswarm plot',
            anchor = 'comp_qm_mulliken_beeswarm',
            plot = beeswarm.plot(self.comp_qm_data['mulliken'], headers, config)
        )

    def comp_qm_geometry_descriptor_chart (self,geometry_descriptor="bonds",min=0.0,max=2.0):
        """ Make the geometry_decriptor section tables and plots"""
        headers = OrderedDict()
        sample_keys=list(self.comp_qm_data[geometry_descriptor].keys())

        for key in self.comp_qm_data[geometry_descriptor][sample_keys[0]].keys(): # hack to get keys from first data set
            headers[key] = {
                'title': key,
                'description': 'Estimated Mulliken Charge',
                'suffix': '',
                'scale': 'Spectral', # 'RdBu',
                'dmin': min,
                'dmax': max,
                'ceiling': max,
                'floor': min,  # known issue with a negative range
                'format':'{:,.2f}',
                'shared_key':'torsion_range'
            }

        log.debug(str(headers))

        # Config for the plot
        config = {
            'namespace': 'Geometry',
            'id': 'comp_qm_geometry',
            'title': 'Geometry'
        }
        self.add_section (
            name = geometry_descriptor,
            anchor = 'comp_qm_geometry',
            plot = table.plot(self.comp_qm_data[geometry_descriptor], headers, config)
        )
        #.. using dmin to ensure correct colour scales for negative ranges, note that min is still required for the beeswarm plot and that behaviour cancels out (?)
        #.. TLDR, the benefit of specifying a range is to keep all the graphs consistent. Maybe require an automated range check, rather than hardcoded. Try without if you like.
        for key in self.comp_qm_data[geometry_descriptor][sample_keys[0]].keys(): # hack to get keys from first data set
            headers[key] = {
                'title': key,
                'description': 'Estimated Mulliken Charge',
                'suffix': '',
                'scale': 'Spectral', # 'RdBu',
                'min': min,
                'max': max,
                'ceiling': max,
                'floor': min,  # known issue with a negative range
                'format':'{:,.2f}',
                'shared_key':'torsion_range'
            }
        self.add_section (
            name = geometry_descriptor,
            anchor = 'comp_qm_geometry',
            plot = beeswarm.plot(self.comp_qm_data[geometry_descriptor], headers, config)
        )
