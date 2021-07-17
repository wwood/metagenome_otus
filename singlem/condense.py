import logging
import tempfile
import csv
import numpy as np
import extern

from queue import PriorityQueue
from collections import OrderedDict
from .singlem_package import SingleMPackage

class Condenser:
    """ Combines otu table output for each marker into a single otu table"""
    @staticmethod
    def condense(**kwargs):
        input_otu_table = kwargs.pop('input_otu_table')
        singlem_packages = kwargs['singlem_packages']
        trim_percent = kwargs.pop('trim_percent') / 100
        output_otu_table = kwargs.pop('output_otu_table')
        krona_output_file = kwargs.pop('krona')
        
        markers = {} # list of markers used
        marker_set = set() # set of markers to check if all markers are included
        samples = {} # otu tables may contain more than one sample ~ {sampleID:{gene:wordtree}}}
        taxon_counter = {} # count number of genes for each taxonomy before removal
        queues = {} # PriorityQueues for processing taxonomy levels in ascending order
        target_domains = {"Archaea": 0, "Bacteria": 0, "Eukaryota": 0}
        
        count = 0
        for path in singlem_packages:
            count += 1
            spkg = SingleMPackage.acquire(path)
            logging.debug("Loading SingleM package: {}".format(spkg.graftm_package_basename()))
            # ensure v3 packages
            if not spkg.version in [3]:
                raise Exception("Only works with v3 packages.")
            markers[spkg.graftm_package_basename()] = spkg.target_domains()
            # count number of markers for each domain
            for domain in spkg.target_domains():
                if domain == "Archaea":
                    target_domains["Archaea"] += 1
                elif domain == "Bacteria":
                    target_domains["Bacteria"] += 1
                elif domain == "Eukaryota":
                    target_domains["Eukaryota"] += 1
                else:
                    raise Exception("Domain: {} not supported.".format(domain))
        logging.info("Loaded %i SingleM packages." % count)
                
        for domain in target_domains:
            if target_domains[domain] in [1, 2]:
                raise Exception("Number of markers for all domains must either be >= 3 or equal to 0. Only {} markers for domain '{}' found".format(target_domains[domain], domain))
        
        excluded_markers = set()
        
        logging.info("Reading input OTU table.")
        with open(input_otu_table) as csvfile:
            reader = csv.DictReader(csvfile, delimiter = "\t")
            for line in reader:
                gene = line["gene"]
                sample = line["sample"]
                taxonomy = line["taxonomy"]
                coverage = float(line["coverage"])
                if gene not in markers:
                    if gene not in excluded_markers:
                        logging.info("Gene: {} not in SingleM packages, excluding hits from this marker...".format(gene))
                        excluded_markers.add(gene)
                    continue
                tax_split = taxonomy.split('; ')
                # ensure OTU is assigned to the domain level or higher
                if len(tax_split) < 2:  # contains at least Root; d__DOMAIN
                    continue
                domain = tax_split[1].strip('d__')
                # ensure the domain of this OTU is targeted by the gene
                if domain not in markers[gene]: 
                    continue
                if sample not in samples:
                    samples[sample] = {}
                    taxon_counter[sample] = {}
                    queues[sample] = PriorityQueue()
                if gene not in marker_set:
                    marker_set.add(gene)
                if gene not in samples[sample]:
                    # create new tree for this marker
                    samples[sample][gene] = WordNode(None, "Root")
                # this taxonomy has been added already?
                taxbuilder = "Root"
                priority = 0
                # process taxonomy levels in ascending order
                for tax in tax_split[1:]:
                    taxbuilder += "; " + tax
                    priority += 1
                    if taxbuilder not in taxon_counter[sample]:
                        queues[sample].put((priority, taxbuilder))
                        taxon_counter[sample][taxbuilder] = set()
                # add to tree
                samples[sample][gene].add_words(tax_split, coverage)
                taxon_counter[sample][taxonomy].add(gene)
        
        for gene in markers:
            if gene not in marker_set:
                logging.warning("Marker %s not found in OTU table, removing..." % gene)
                for domain in markers[gene]:
                    target_domains[domain] -= 1
                    if target_domains[domain] in [1, 2]:
                        raise Exception("Number of markers for all domains must either be >= 3 or equal to 0. Only {} markers for domain '{}' found".format(target_domains[domain], domain))
        
        ### Trim the taxonomy coverages
        logging.info("Trimming constructed taxonomy trees.")
        for sample in queues:
            while not queues[sample].empty():
                taxonomy = queues[sample].get()[1] # tuple: (priority, taxonomy)
                tax_split = taxonomy.split("; ")
                domain = tax_split[1].strip("d__")
                genes = {}
                for gene in samples[sample]:
                    node = samples[sample][gene].get_tree(taxonomy.split("; "))
                    # node returns False if it does not exist
                    if node:
                        cov = node.get_full_coverage()
                        genes[gene] = cov
                # do not trim if taxonomy only remains in 2 markers or less
                if len(genes) > 2:
                    trimmed_mean = _tmean(list(genes.values()), trim_percent)
                    # construct list of differences (diff, marker)
                    diffs = [(abs(trimmed_mean-genes[gene]), gene) for gene in genes]
                    diff_tail = _tail(diffs, trim_percent)
                    for diff in diff_tail:
                        # remove taxonomy from markers with largest tails
                        gene = diff[1]
                        samples[sample][gene].remove_tree(tax_split)
                        logging.debug("Removed taxonomy node %s for gene %s." % (taxonomy, gene))
                else:
                    # push coverage down since we cannot be confident at this level
                    i = 1
                    while i < len(tax_split)-1:
                        # find taxonomy that has at least 3 hits across markers, otherwise place at domain level
                        new_tax = "; ".join(tax_split[:-i])
                        if len(taxon_counter[sample][new_tax]) >= 3:
                            break
                        i += 1
                    for gene in genes:
                        node = samples[sample][gene].get_tree(tax_split)
                        if node.is_key():
                            samples[sample][gene].get_tree(new_tax.split("; ")).coverage += node.coverage
                            node.coverage = 0
                            taxon_counter[sample][new_tax].add(gene)
                            taxon_counter[sample][taxonomy].remove(gene)
                    
        ### Condense trees into final OTU tables
        logging.info("Condensing trees into OTU tables")
        fieldnames = ["sample", "coverage", "taxonomy"]
        sample_tempfiles = []
        cmd = 'ktImportText -o %s' % krona_output_file
        
        with open(output_otu_table, 'w+') as final_otu:
            csv_writer = csv.DictWriter(final_otu, fieldnames, delimiter = '\t')
            csv_writer.writeheader()
            for sample in taxon_counter:
                # recreate queue
                queueList = sorted([taxonomy for taxonomy in taxon_counter[sample]])
                f = tempfile.NamedTemporaryFile(prefix='singlem_condense_krona', mode= 'w')
                sample_tempfiles.append(f)
                for taxonomy in queueList:
                    domain = taxonomy.split("; ")[1].strip("d__")
                    # ensure this taxonomy has been found in at least 3 markers
                    if len(taxon_counter[sample][taxonomy]) < 3:
                        continue
                    count = 0
                    otu = OrderedDict({"sample": sample, "coverage": 0, "taxonomy": taxonomy})
                    for marker in taxon_counter[sample][taxonomy]:
                        node = samples[sample][marker].get_tree(taxonomy.split("; "))
                        if node:
                            if node.is_key():
                                otu["coverage"] += node.coverage
                                count += 1
                    # calculate otu coverage
                    if otu["coverage"] == 0 or count == 0:
                        continue
                    otu["coverage"] = str((otu["coverage"]) / target_domains[domain])
                    # write otu to otu tables
                    csv_writer.writerow(otu)
                    tax_split = otu["taxonomy"].split('; ')
                    if tax_split[0] == "Root": tax_split = tax_split[1:]
                    f.write('\t'.join([otu["coverage"]] + tax_split))
                    f.write('\n')
                f.flush()
                cmd += " %s,'%s'" % (f.name, sample)
        logging.info("Writing krona %s" % krona_output_file)
        extern.run(cmd)
        for f in sample_tempfiles:
            f.close()
        logging.info("Finished")

def _tmean(data, proportiontocut):
    """
    Returns trimmed mean of an array. 
    
    Slices at least the first and last values of the array.
    """
    a = sorted(data)
    cut = int(np.floor(len(a) * proportiontocut))
    if cut == 0:
        cut = 1
    return np.mean(a[cut:-cut])

def _tail(data, proportiontocut):
    """
    Returns the tail slice of an array. 
    
    Slice will contain at least the final 2 values.
    """
    a = sorted(data)
    cut = int(np.floor(len(a) * proportiontocut))
    if cut < 2:
        cut = 2
    return a[-cut:]

class WordNode:
    def __init__(self, parent, word):
        self.parent = parent # WordNode object
        self.word = word
        self.children = {}
        self.coverage = 0
        
    def is_key(self):
        if self.coverage == 0:
            return False
        return True
    
    def add_words(self, word_list, coverage):
        """
        Add coverage to this taxonomy to tree.

        Parameters
        ----------
        word_list : LIST
            list of taxonomy terms to be parsed.
        coverage : FLOAT
            coverage to be inserted at the tail end of word_list.

        Raises
        ------
        Exception
            if the current word in word_list mismatches the current node.
        """
        # first word is current
        if self.word == word_list[0]:
            if len(word_list) > 1:
                if word_list[1] in self.children:
                    self.children[word_list[1]].add_words(word_list[1:], coverage)
                else:
                    # word is not present, add to children
                    self.children[word_list[1]] = WordNode(self, word_list[1])
                    self.children[word_list[1]].add_words(word_list[1:], coverage)
            else:
                # last word
                self.coverage += coverage
        else:
            raise Exception("Word %s does not match current node %s" % (word_list[0]), self.word)
    
    def get_tree(self, word_list):
        """ 
        Gets node associated with this taxonomy. 
        """
        if word_list[0] == self.word:
            if len(word_list) > 1:
                if word_list[1] in self.children:
                    return self.children[word_list[1]].get_tree(word_list[1:])
            else:
                return self
        # tree not present
        return False
    
    def remove_tree(self, word_list):
        """ 
        Removes node and descendants described by the word list
        from this tree. 
        """
        parent_node = self.get_tree(word_list[:-1])
        child = word_list[-1]
        if parent_node:
            if child in parent_node.children:
                # delete this child
                parent_node.children.pop(child)
    
    def get_full_coverage(self):
        """
        Get full coverage sum of this node and all its children.
        """
        covs = 0
        if self.is_key():
            covs += self.coverage
        if len(self.children) == 0: # leaf node
            return covs
        for node in self.children.values():
            covs += node.get_full_coverage()
        return covs
