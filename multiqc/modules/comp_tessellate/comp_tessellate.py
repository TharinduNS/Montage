#!/usr/bin/env python

""" MultiQC module to parse output from Tessellate """

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
    """Tessellate module """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Tessellate', anchor='comp_tessellate',
        href="http://bitbucket.org/scientificomputing.com",
        info="A program for tesselating ring cycles, identifying conformations from coordinates and molecular trajectories")

        # Set up data structure
        #self.comp_tessellate_data = dict()
        self.comp_tessellate_data = {
            'five': {},
            'six': {},
            'seven': {},
            'eight': {},
            'macro': {},
            'five_numeric': {},#. TODO BUGFIX data duplication and code whiff to be removed 
            'six_numeric': {},#. TODO BUGFIX data duplication and code whiff to be removed 
            'seven_numeric': {},#. TODO BUGFIX data duplication and code whiff to be removed 
            'eight_numeric': {},#. TODO BUGFIX data duplication and code whiff to be removed
            'all': {}, #. TODO BUGFIX data duplication and code whiff to be removed
            'all_json': {} #TODO BUGFIX data duplication and code whiff to be removed 
        }


        #. decide whether to create subsamples based on pdbid (bad idea for timeseries, good idea for PDB structures)
        self.subsamples=False

        for f in self.find_log_files('comp_tessellate', filehandles=True):
            self.parse_comp_tessellate_log(f)

        # Filter to strip out ignored sample names
        self.comp_tessellate_data = self.ignore_samples(self.comp_tessellate_data)

        if len(self.comp_tessellate_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.comp_tessellate_data)))

        # Write parsed report data to a file
        self.write_data_file(self.comp_tessellate_data, 'multiqc_comp_tessellate')

        # Simple Time Series Plot
        if len(self.comp_tessellate_data['all'])>0:
            self.add_section( plot = self.comp_tessellate_timeseries_plot() )
        if len(self.comp_tessellate_data['five'])>0:
            self.comp_tessellate_conformer_chart(ring='five')
            self.comp_tessellate_conformer_scatter(ring='five')
        if len(self.comp_tessellate_data['six'])>0:
            self.comp_tessellate_conformer_chart(ring='six')
            self.comp_tessellate_conformer_scatter(ring='six')
        if len(self.comp_tessellate_data['seven'])>0:
            self.comp_tessellate_conformer_chart(ring='seven')
            self.comp_tessellate_conformer_scatter(ring='seven')
        if len(self.comp_tessellate_data['eight'])>0:
            self.comp_tessellate_conformer_chart(ring='eight')
            self.comp_tessellate_conformer_scatter(ring='eight')
        if len(self.comp_tessellate_data['macro'])>0:
            self.comp_tessellate_conformer_chart(ring='macro')

        #if len(self.comp_tessellate_data['all_json'])>0:
        #    self.comp_tessellate_conformer_table()


    def parse_comp_tessellate_log(self, f):
        s_name = None
        conformers={}
        #numerics_all_size=[]
        numerics_all_sizes={'5':[],'6':[],'7':[],'8':[],'macro':[]} #duplication to solve quickly, TODO BUGFIX
        conformers_all_sizes={'5':{},'6':{},'7':{},'8':{},'macro':{}} #duplication to solve quickly, TODO BUGFIX
        pdb_all_sizes={'5':{},'6':{},'7':{},'8':{},'macro':{}} #duplication to solve quickly, TODO BUGFIX

        s_name = self.clean_s_name(f['s_name'],f['root'])
        #. json file example
        import json
        firstline = f['f'].readline().split() # check the format
        if len(firstline)==3:
            _,tessellate_version,fileformat = firstline
        elif len(firstline)==4:
            _,tessellate_version,fileformat,input_format = firstline
            if "pdb" in input_format:
                self.subsamples=True
        else:
            log.debug("Unable to read firstline of file")
        log.debug("Version %s %s", tessellate_version,fileformat)
        if fileformat=="txt":
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
        elif fileformat=="json":
            all_json_data=json.loads("".join(f['f'].readlines()))
            log.debug(all_json_data)
            #. loop over json list
            for itm in all_json_data:
                conformer=itm["conformer"]
                ringsize=itm["ringsize"]
                numeric=itm["numeric"]
                pdbid=itm["pdbid"]
                try:
                    #. add data to numerics no matter what
                    #.. Add a colorbrewer2.org scheme. 5-class Set3 for qualitative
                    if 'E' in conformer:
                        ringcolor="#8dd3c7"
                    elif 'T' in conformer:
                        ringcolor="#ffffb3"
                    elif 'UAP' in conformer:
                        ringcolor="#bebada"
                    elif 'P' in conformer:
                        ringcolor="#fb8072"
                    else:
                        ringcolor="#80b1d3"
                    numerics_all_sizes[str(ringsize)].append({'x':len(numerics_all_sizes[str(ringsize)]),'y':numeric,'name':conformer, 'color':ringcolor})
                    if conformer in conformers.keys():
                        log.debug("Add to conformer %s %s",conformer,conformers.keys())
                        conformers[conformer]+=1
            
                    else: #initialise count to 1 for this conformer
                        log.debug("Initialise conformer %s %s",conformer,conformers.keys())
                        conformers[conformer]=1
                    if conformer in conformers_all_sizes[str(ringsize)].keys():
                        conformers_all_sizes[str(ringsize)][conformer]+=1
            
                    else: #initialise count to 1 for this conformer
                        conformers_all_sizes[str(ringsize)][conformer]=1

                    #. stratify by pdbid and start wondering if something like pandas could do this better...
                    if self.subsamples:
                        if pdbid in pdb_all_sizes[str(ringsize)].keys():
                            if conformer in pdb_all_sizes[str(ringsize)][pdbid].keys():
                                pdb_all_sizes[str(ringsize)][pdbid][conformer]+=1
                            else:
                                pdb_all_sizes[str(ringsize)][pdbid][conformer]=1
                        else:
                            pdb_all_sizes[str(ringsize)][pdbid]={conformer:1}

                except Exception as e:
                    log.debug("Except conformer not included in dictionary %s", conformer)
                    raise e
        else:
            log.error("Unknown input file format %s",fileformat)


        #. add data to the tessellate section - all json data
        self.add_data_source(f, section='all_json') #section='mulliken')
        datadict = {str(idx): dict(key) for idx, key in enumerate(all_json_data)}
        self.comp_tessellate_data['all_json'][s_name] = datadict

        #. add data to the tessellate section
        self.add_data_source(f, section='all') #section='mulliken')
        self.comp_tessellate_data['all'][s_name] = conformers
        log.debug("Data added for this log: %s",str(self.comp_tessellate_data['all'][s_name]))
        #. add data to the tessellate 5 section
        self.add_data_source(f, section='five') #section='mulliken')
        if len(conformers_all_sizes['5'])>0:
            self.comp_tessellate_data['five'][s_name] = conformers_all_sizes['5']
            log.debug("Data added for this log: %s",str(self.comp_tessellate_data['five'][s_name]))
        #. add data to the tessellate 6 section
        self.add_data_source(f, section='six') #section='mulliken')
        if len(conformers_all_sizes['6'])>0:
            self.comp_tessellate_data['six'][s_name] = conformers_all_sizes['6']
            log.debug("Data added for this log: %s",str(self.comp_tessellate_data['six'][s_name]))
        #. add data to the tessellate 7 section
        self.add_data_source(f, section='seven') #section='mulliken')
        if len(conformers_all_sizes['7'])>0:
            self.comp_tessellate_data['seven'][s_name] = conformers_all_sizes['7']
            log.debug("Data added for this log: %s",str(self.comp_tessellate_data['seven'][s_name]))
        #. add data to the tessellate 8 section
        self.add_data_source(f, section='eight') #section='mulliken')
        if len(conformers_all_sizes['8'])>0:
            self.comp_tessellate_data['eight'][s_name] = conformers_all_sizes['8']
            log.debug("Data added for this log: %s",str(self.comp_tessellate_data['eight'][s_name]))
        #. add numerics data for scatter plot
        self.comp_tessellate_data['five_numeric'][s_name] = numerics_all_sizes['5']
        self.comp_tessellate_data['six_numeric'][s_name] = numerics_all_sizes['6']
        self.comp_tessellate_data['seven_numeric'][s_name] = numerics_all_sizes['7']
        self.comp_tessellate_data['eight_numeric'][s_name] = numerics_all_sizes['8']
        #. segregate data by pdbid as well....
        #. loop over all pdb data size and each key
        if self.subsamples:
            for pdbkey in pdb_all_sizes['5']:
                self.comp_tessellate_data['five']["_".join([s_name,pdbkey])] = pdb_all_sizes['5'][pdbkey]
            for pdbkey in pdb_all_sizes['6']:
                self.comp_tessellate_data['six']["_".join([s_name,pdbkey])] = pdb_all_sizes['6'][pdbkey]
            for pdbkey in pdb_all_sizes['7']:
                self.comp_tessellate_data['seven']["_".join([s_name,pdbkey])] = pdb_all_sizes['7'][pdbkey]
            for pdbkey in pdb_all_sizes['8']:
                self.comp_tessellate_data['eight']["_".join([s_name,pdbkey])] = pdb_all_sizes['8'][pdbkey]
    


        s_name = None
        conformers = {}

    def comp_tessellate_timeseries_plot (self):
        """ Make the HighCharts HTML to plot the timeseries """

        # Specify the order of the different possible categories
        headers = OrderedDict()
        sample_keys=list(self.comp_tessellate_data['all'].keys())
        log.debug("Sample keys: %s", sample_keys)

        #. For ring pucker the first sample may/may not have all keys. Aggregate keys
        all_pucker_keys_found=[]
        for key in sample_keys:
            all_pucker_keys_found.extend(list(self.comp_tessellate_data['all'][key].keys()))

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
                'shared_key':'tessellate'
            }

        # Config for the plot
        config = {
            'id': 'comp_tessellate_count',
            'title': 'Tessellate: Pucker distribution',
            'ylab': '# Count',
        }
        return bargraph.plot(self.comp_tessellate_data['all'], headers, config)


    def comp_tessellate_conformer_chart (self,ring="five",min=0.0,max=2.0):
        """ Make the bar chart for any ring size"""
        headers = OrderedDict()
        sample_keys=list(self.comp_tessellate_data[ring].keys())
        log.debug("Sample keys: %s", sample_keys)
        #. For ring pucker the first sample may/may not have all keys. Aggregate keys
        all_pucker_keys_found=[]
        for key in sample_keys:
            all_pucker_keys_found.extend(list(self.comp_tessellate_data[ring][key].keys()))
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
                'shared_key':'tessellate'
            }

        # Config for the plot
        created_title =  'Tessellate: Pucker distribution for ' + ring
        config = {
            'id': 'comp_tessellate_count_by_size',
            'title': created_title, #'Tessellate: Count of Puckers by size and by canonical conformer',
            'namespace': 'Tessellate',
            'ylab': '# Count',
        }
        self.add_section (
            name = ring,
            anchor = 'comp_tessellate_count_by_size',
            plot = bargraph.plot(self.comp_tessellate_data[ring], headers, config)
        )

    def comp_tessellate_conformer_scatter (self,ring="five",min=0.0,max=2.0):
        """ Make the scatter chart for any ring size"""
        headers = OrderedDict()
        ring_key=ring+"_numeric"
        sample_keys=list(self.comp_tessellate_data[ring_key].keys())
        log.debug("Sample keys: %s", sample_keys)
        log.debug("Sample keys: %s", self.comp_tessellate_data[ring_key])
        #. For ring pucker the first sample may/may not have all keys. Aggregate keys
        # Config for the plot
        created_title =  'Tessellate: Pucker series for ' + ring
        config = {
            'id': 'comp_tessellate_count_by_size',
            'title': created_title, 
            'namespace': 'Tessellate',
            'xlab': 'Time or Sequence order',
            'ylab': 'Numeric Pucker ID',
        }
        created_name =  ring+' series'
        self.add_section (
            name = created_name,
            anchor = 'comp_tessellate_count_by_size',
            plot = scatter.plot(self.comp_tessellate_data[ring_key], config)
        )

    def comp_tessellate_conformer_table (self,min=0.0,max=2.0):
        """ Make a table of all data """
        headers = OrderedDict()
        sample_keys=list(self.comp_tessellate_data['all_json'].keys())
        log.debug("Sample keys: %s", sample_keys)
        log.critical("all_json %s", self.comp_tessellate_data['all_json']) 

        #make an assumption that all json elements have the same headings
        for key in self.comp_tessellate_data['all_json'][sample_keys[0]].keys(): # hack to get keys from first data set
            headers[key] = {
                'title': key,
                'description': 'Ring pucker data',
                'suffix': '',
                'scale': 'Spectral', # 'RdBu',
                'dmin': min,
                'dmax': max,
                'ceiling': max,
                'floor': min,  # known issue with a negative range
                'format':'{:,.2f}',
                'shared_key':'tessellate'
            }

        # Config for the plot
        created_title =  'Tessellate: Pucker table'
        config = {
            'id': 'comp_tessellate_table',
            'title': created_title, #'Tessellate: Count of Puckers by size and by canonical conformer',
            'namespace': 'Tessellate',
        }
        self.add_section (
            name = 'Pucker table',
            anchor = 'comp_tessellate_table',
            plot = table.plot(self.comp_tessellate_data['all_json'], headers, config)
        )

