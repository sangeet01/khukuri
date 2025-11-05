"""Fragment library for drug design"""

import logging

logger = logging.getLogger('khukuri')


class FragmentLibrary:
    """Fragment-based drug design library"""
    
    def __init__(self):
        self.scaffolds = {
            'benzene': 'c1ccccc1',
            'pyridine': 'c1ccncc1',
            'furan': 'c1ccoc1',
            'thiophene': 'c1ccsc1',
            'pyrrole': 'c1cc[nH]c1',
            'imidazole': 'c1cnc[nH]1'
        }
        
        self.substituents = {
            'halogens': ['F', 'Cl', 'Br'],
            'alkyl': ['C', 'CC', 'CCC'],
            'amino': ['N', 'NC'],
            'hydroxyl': ['O'],
            'carbonyl': ['C(=O)C']
        }
    
    def get_scaffolds(self, scaffold_type='all'):
        """Get scaffold SMILES"""
        if scaffold_type == 'all':
            return list(self.scaffolds.values())
        return [self.scaffolds.get(scaffold_type, 'c1ccccc1')]
    
    def get_substituents(self, sub_type='all'):
        """Get substituent SMILES"""
        if sub_type == 'all':
            return [s for subs in self.substituents.values() for s in subs]
        return self.substituents.get(sub_type, ['C'])
