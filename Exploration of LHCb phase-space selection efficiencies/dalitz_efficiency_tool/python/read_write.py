from ROOT import RDF, TGenPhaseSpace
import numpy as np
import uproot
from python.helper_funcs_polynomial import find_pairs

def generate_flat_event(generator):
    """
    generates flat events for a given instance of TGenPhaseSpace

    :param generator: configured instance of TGenPhaseSpace:Generator
    :returns: flat event generator instance
    """
    while True:
        weight = generator.Generate()
        height = np.random.uniform(0, 1.0/weight)
        if height < 1.0:
            return generator

def write_file(event_data_obj, file_name, tree_name='D0_rest'):
    rdf = RDF.MakeNumpyDataFrame({
        'm12': event_data_obj.m12_array,
        'm13': event_data_obj.m13_array
    })

    rdf.Snapshot(tree_name, f'event_data/{file_name}.root')

def read_file(file_name, tree_name='D0_rest'):
    with uproot.open(f'event_data/{file_name}.root') as f:
        tree = f[tree_name]
        m12_array = np.array(tree['m12'].array())
        m13_array = np.array(tree['m13'].array())

    event_data = Events(m12_array, m13_array)

    return event_data

class Events:

    def __init__(self, m12_array=None, m13_array=None):
        self.m12_array = m12_array
        self.m13_array = m13_array

    def generate_events(self, parent_vect, masses, no_decays):

        event = TGenPhaseSpace()                                        # Uses the same cycling seed every run
        event.SetDecay(parent_vect, 3, masses)

        self.m12_array = np.array([])
        self.m13_array = np.array([])

        count = 0
        print()

        for _ in range(no_decays):
            event = generate_flat_event(event)

            m1_vect = event.GetDecay(0)

            m2_vect = event.GetDecay(1)
            m3_vect = event.GetDecay(2)

            m12_vect = m1_vect + m2_vect
            m13_vect = m1_vect + m3_vect

            self.m12_array = np.append(self.m12_array, m12_vect.M2())
            self.m13_array = np.append(self.m13_array, m13_vect.M2())

            count += 1
            if count % (no_decays / 10) == 0:
                print(f"Generated {count} events...")

    def get_arrays(self):

        return self.m12_array, self.m13_array

    def artificial_rejection_uniform(self):

        rejected_self = self

        i = np.random.uniform(size=len(self.m13_array))
        j = np.random.uniform(size=len(self.m13_array))
        events_remove_ind = np.where(i > j)[0]

        rejected_self.m13_array = np.delete(self.m13_array, events_remove_ind)
        rejected_self.m12_array = np.delete(self.m12_array, events_remove_ind)

        return rejected_self

    def artificial_rejection_lin(self):

        rejected_self = self

        # This assumes that if we for example fail to 'detect' m3, we ignore the entire data entry even if we have m1 & m2

        r = np.random.uniform(np.min(self.m13_array), np.max(self.m13_array), size=len(self.m13_array))
        events_remove_ind = np.where(r * 1.25 < self.m13_array)[0]

        rejected_self.m13_array = np.delete(self.m13_array, events_remove_ind)
        rejected_self.m12_array = np.delete(self.m12_array, events_remove_ind)

        return rejected_self

    def artificial_rejection_quad(self):

        rejected_self = self

        m13_sqrt = np.sqrt(self.m13_array)
        r = np.random.uniform(np.min(m13_sqrt), np.max(m13_sqrt), size=len(m13_sqrt))
        events_remove_ind = np.where(r < m13_sqrt)[0]

        rejected_self.m13_array = np.delete(self.m13_array, events_remove_ind)
        rejected_self.m12_array = np.delete(self.m12_array, events_remove_ind)

        return rejected_self

    def artificial_rejection_cubic2d(self, params):

        # params should be according to the find_paris ordering used in the MC fits

        # ordering for degree==3 is : [[0, 0], [1, 0], [2, 0], [3, 0], [0, 1], [1, 1], [2, 1], [0, 2], [1, 2], [0, 3]]

        pairs = find_pairs(3)
        temp = np.array([np.power(self.m12_array, pair[0]) * np.power(self.m13_array, pair[1]) for pair in pairs])
        f_i = np.dot(params, temp)
        f_max = np.max(f_i)
        r_i = np.random.uniform(0, 1, size=len(f_i))

        events_remove_ind = np.where(r_i > f_i/f_max)[0]
        redacted_m12_array = np.delete(self.m12_array, events_remove_ind)
        redacted_m13_array = np.delete(self.m13_array, events_remove_ind)

        return Events(redacted_m12_array, redacted_m13_array)
