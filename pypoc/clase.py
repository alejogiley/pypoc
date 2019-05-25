# point of this class
class _Acyls:
        """Comments"""

        def __init__(self, resid, bondd, angle):
                # resid: residue id
                # bonds: bond id id
                # resid: residue id
                # angle: angle in rad
                self.res = [resid]
                self.bon = [bondd]
                self.ang = [angle]
                # order
                self.ods = self.order2
                # counter
                self.num = 1

        def add_data(self, resid, bondd, angle):
                # increment 
                self.num += 1
                # fill list
                self.res.append(resid)
                self.bon.append(bondd)
                self.ang.append(angle)

        @property
        def order2(self):
                """Comment"""
                # set function
                f = lambda x: 0.25 + 0.75 * math.cos(2.0 * x)
                # apply to all angles
                self.ods  = sum(map(f, self.ang))
                # average order
                self.ods /= len(self.ang)
                return self.ods

        @staticmethod
        def average(x):
                """Comment"""
                # average angle sine
                ssin = sum(map(math.sin, x))
                # average angle cosine
                scos = sum(map(math.cos, x))
                # average angle
                return math.atan2(ssin, scos)

# point of this function
def add_to_dict(data, ids, resids, bonds, angles):
        """ Checks wether an id has been added to data dict.
            If it has been added, increment record,
            otherwise, creates a new record. """

        for c,v in enumerate(ids):
                if data.get(v, None) is None:
                        # create object
                        data[v] = _Acyls(resids[c], bonds[c], angles[c])
                else:
                        # append to object
                        data[v].add_data(resids[c], bonds[c], angles[c])


