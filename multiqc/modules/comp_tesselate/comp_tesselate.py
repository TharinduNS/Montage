#!/usr/bin/env python

""" MultiQC module to parse output from Tesselate """

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
    """Tesselate module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Tesselate', anchor='comp_tesselate',
        href="http://bitbucket.org/scientificomputing.com",
        info="A program for tesselating ring cycles, identifying conformations from coordinates and molecular trajectories")

        # Set up data structure
        self.comp_tesselate_data = dict()

        for f in self.find_log_files('comp_tesselate', filehandles=True):
            self.parse_comp_tesselate_log(f)

        # Filter to strip out ignored sample names
        self.comp_tesselate_data = self.ignore_samples(self.comp_tesselate_data)

        if len(self.comp_tesselate_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.comp_tesselate_data)))

        # Write parsed report data to a file
        self.write_data_file(self.comp_tesselate_data, 'multiqc_comp_tesselate')

        # Simple Time Series Plot
        self.add_section( plot = self.comp_tesselate_timeseries_plot() )

    def parse_comp_tesselate_log(self, f):
        s_name = None
        conformers={}
        s_name = self.clean_s_name(f['s_name'],f['root'])
        if s_name is not None:
            for l in f['f']:
                chunked = l.split(" ")
                log.debug(chunked)
                if len(chunked)>1:
                    #TODO improve use a sed or regexp
                    conformer= chunked[1].strip("'").strip("(").strip("'").strip(",").strip("'")
                    log.debug(conformer)
                    try:
                        if conformer in conformers.keys():
                            log.debug("Add to conformer %s %s",conformer,conformers.keys())
                            conformers[conformer]+=1
                        else: #initialise count to 1 for this conformer
                            log.debug("Initialise conformer %s %s",conformer,conformers.keys())
                            conformers[conformer]=1
                    except Exception as e:
                        log.debug("Except conformer not included in dictionary %s", conformer)
                        raise e

        #. add data to the tesselate section
        self.add_data_source(f, s_name) #section='mulliken')
        self.comp_tesselate_data[s_name] = conformers
        log.debug("Data added for this log: %s",str(self.comp_tesselate_data[s_name]))

        s_name = None
        conformers = {}

    def comp_tesselate_timeseries_plot (self):
        """ Make the HighCharts HTML to plot the timeseries """

        # Specify the order of the different possible categories
        headers = OrderedDict()
        sample_keys=list(self.comp_tesselate_data.keys())
        log.debug("Sample keys: %s", sample_keys)

        #. For ring pucker the first sample may/may not have all keys. Aggregate keys
        all_pucker_keys_found=[]
        for key in sample_keys:
            all_pucker_keys_found.extend(list(self.comp_tesselate_data[key].keys()))

        log.debug("All pucker keys: %s", all_pucker_keys_found)

        for key in all_pucker_keys_found: # get all_keys
            headers[key] = {
                'title': key,
                'description': 'Ring pucker conformation',
                'suffix': '',
                'scale': 'Spectral', # 'RdBu',
                'dmin': min,
                'dmax': max,
                'ceiling': max,
                'floor': min,  # known issue with a negative range
                'format':'{:,.2f}',
                'shared_key':'tesselate'
            }

        # Config for the plot
        config = {
            'id': 'comp_tesselate_count',
            'title': 'Tesselate: Count of Puckers by canonical conformer',
            'ylab': '# Count',
        }
        return bargraph.plot(self.comp_tesselate_data, headers, config)
