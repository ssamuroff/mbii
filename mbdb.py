import MySQLdb as mdb
import numpy as np
from numpy.core.records import fromarrays
import pylab as plt
import treecorr
plt.switch_backend('pdf')
#plt.style.use('y1a1')


class mbdb(object):

    def __init__(self, sqlserver='localhost', user='flanusse', password='mysqlpass',
                 dbname='mb2_hydro', unix_socket='/home/rmandelb.proj/flanusse/mysql/mysql.sock'):
        try:
            self.db = mdb.connect(sqlserver, user, password, dbname, unix_socket=unix_socket)
        except:
            print("Error when connecting to the database")

    def get_columns(self, table):
        """
        Returns the list of columns in the specificied table
        """
        c = self.db.cursor()
        c.execute('SELECT * FROM %s LIMIT 1;'%table)
        return [i[0] for i in c.description]

    def get(self, table, fields, cond="", fatal_errors=True):
        """
        Returns the sql query as a nice numpy recarray
        expects the list of fields in the format fields='a,b,c'
        """

        sql = "SELECT %s FROM %s WHERE %s;"%(fields, table, cond)
        print sql
        try:
            # prepare a cursor for the query
            cursor = self.db.cursor()
            cursor.execute(sql)
            print("Fetching %d entries" % cursor.rowcount)
        except:
            if fatal_errors:
                    raise
            else:
                print("Error when runnning the SQL command")
                return

        results = fromarrays(np.array(cursor.fetchall()).squeeze().T,names=fields)
        return results

    def cross_match(self, source, table, fields, match_column, match_column2='', fatal_errors=True):

        # Build the SQL query
        if match_column2=='':
            match_column2=match_column
        print 'Will cross match column %s in the table provided with %s in table %s'%(match_column,match_column2,table)
        print 'Building query...'
        sql = "SELECT %s FROM %s WHERE %s IN ("%(fields, table, match_column)
        for row in source[match_column]:
            sql+="'%s',"%str(row)
        sql = sql[:-1] + ')'
        try:
            # prepare a cursor for the query
            cursor = self.db.cursor()
            cursor.execute(sql)
            print("Fetching %d entries" % cursor.rowcount)
        except:
            if fatal_errors:
                    raise
            else:
                print("Error when runnning the SQL command")
                return

        results = fromarrays(np.array(cursor.fetchall()).squeeze().T,names=fields)
        return results

    def get_sql(self, sql, fields):
        try:
            # prepare a cursor for the query
            cursor = self.db.cursor()
            cursor.execute(sql)
            print("Fetching %d entries" % cursor.rowcount)
        except:
            print("Error when runnning the SQL command")
            return
        results = fromarrays(np.array(cursor.fetchall()).squeeze().T,names=fields)
        return results

    def change_snapshot(self,snapshot):
        print "Jumping to snapshot number %d"%snapshot
        self.snapshot = snapshot
        return None



class groups(mbdb):
    def __init__(self, snapshot=85):
        super(groups, self).__init__()

        print "Will load group data for snapshot %d"%(snapshot)
        self.snapshot=snapshot

        #info = self.get(table='subfind_halos' , fields='groupId', cond='subfind_halos.snapnum = %d'%(self.snapshot))
        #self.groups = np.unique(info['groupId'])
        return None

    def get_tomographic_xigg(self, i, j, nbins=8):
        # Define some bins in Rpar
        edges = np.linspace(-2e3,2e3,nbins)
        Rpar = (edges[1:]+edges[:-1])/2

        # Load the data. Again.
        data1 =  self.get(table='subfind_halos' , fields='groupId, x, y, z', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, i ))
        data2 =  self.get(table='subfind_halos' , fields='groupId, x, y, z', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, j))

        vec = []
        rvec=[]
        # Calculate xi(Rperp) in each bin of Rpar
        for j, (lower,upper) in enumerate(zip(edges[:-1], edges[1:])):
            print j, lower, upper
            Rper,w,werr = self.calc_xi_perp(data1, data2, min_rpar=lower, max_rpar=upper)
            vec.append(w)
            rvec.append(Rper)

        return Rpar, rvec, vec

#    def make_contours(self, R, xi_grid):
#        Pi = 
#        xx,yy = np.meshgrid()
#C = plt.contour(Pi,xx.T,res, 8, colors='purple',linestyles='-', linewidth=.5)
#



    def calc_xi_perp(self, data1, data2,  min_rpar, max_rpar, nbins=20, slop=0.1, randoms=True):


        # Build a catalogue of random points drawn from the same volume
        rx = np.random.random(size=data1['x'].size) * (data1['x'].max()-data1['x'].min()) + data1['x'].mean()
        ry = np.random.random(size=data1['x'].size) * (data1['y'].max()-data1['y'].min()) + data1['y'].mean()
        rz = np.random.random(size=data1['x'].size) * (data1['z'].max()-data1['z'].min()) + data1['z'].mean()

        # Create the catalogues
        cat_i = treecorr.Catalog(x=data1['x'], y=data1['y'], z=data1['z'])
        cat_j = treecorr.Catalog(x=data2['x'], y=data2['y'], z=data2['z'])
        rancat_i  = treecorr.Catalog(x=rx, y=ry, z=rz)
        rancat_j  = treecorr.Catalog(x=rx, y=ry, z=rz)

        nn = treecorr.NNCorrelation(nbins=nbins, min_rpar=min_rpar, max_rpar=max_rpar, min_sep=15, max_sep=10e3, bin_slop=slop)
        rn = treecorr.NNCorrelation(nbins=nbins, min_rpar=min_rpar, max_rpar=max_rpar, min_sep=15, max_sep=10e3, bin_slop=slop)
        nr = treecorr.NNCorrelation(nbins=nbins, min_rpar=min_rpar, max_rpar=max_rpar, min_sep=15, max_sep=10e3, bin_slop=slop)
        rr = treecorr.NNCorrelation(nbins=nbins, min_rpar=min_rpar, max_rpar=max_rpar, min_sep=15, max_sep=10e3, bin_slop=slop)
        nn.process(cat_i,cat_j, metric='Rperp') #, metric='Periodic')
        rn.process(rancat_i,cat_j, metric='Rperp') #, metric='Periodic')
        nr.process(cat_i,rancat_j, metric='Rperp') #, metric='Periodic')
        rr.process(rancat_i,rancat_j, metric='Rperp') #, metric='Periodic')

        R = np.exp(nn.meanlogr)
        if randoms:
            w, werr = nn.calculateXi(rr,dr=nr,rd=rn)
        else:
            w, werr = nn.calculateXi(rr,dr=None,rd=None)
        werr = np.sqrt(werr)

        return R, w, werr

    def calc_gg(self, data1, data2, type1='central', type2='satellite', nbins=20, min_sep=44, max_sep=6e3, slop=0.1, weights=None):
        fn = getattr(self, '_calc_gg_%s_%s'%(type1,type2))
        print 'Will generate %s %s correlation'%(type1,type2)
        return fn(data1, data2, nbins=20, min_sep=44, max_sep=6e3, slop=0.1, weights=None)

    def _calc_gg_central_central(self, data1, data2, nbins=20, min_sep=44, max_sep=6e3, slop=0.1, randoms=True, weights=None, return_all=False):
        
        mask1 = (data1['central']==1)
        mask2 = (data2['central']==1)

        # Build a catalogue of random points drawn from the same volume
        rx_j = np.random.random(size=data2['x'][mask2].size) * (data2['x'][mask2].max()-data1['x'][mask1].min()) + data2['x'][mask2].mean()
        ry_j = np.random.random(size=data2['x'][mask2].size) * (data2['y'][mask2].max()-data1['y'][mask1].min()) + data2['y'][mask2].mean()
        rz_j = np.random.random(size=data2['x'][mask2].size) * (data2['z'][mask2].max()-data1['z'][mask1].min()) + data2['z'][mask2].mean()

        # Create the catalogues
        cat_i = treecorr.Catalog(w=weights, x=data1['x'][mask1], y=data1['y'][mask1], z=data1['z'][mask1])
        cat_j = treecorr.Catalog(w=weights, x=data2['x'][mask2], y=data2['y'][mask2], z=data2['z'][mask2])
        rancat_j  = treecorr.Catalog(x=rx_j, y=ry_j, z=rz_j)

        # Compute
        nn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
        nn.process(cat_i,cat_j) #, metric='Periodic')

        # Get the minimum variance estimator
        R = np.exp(nn.meanlogr)
        w = (nn.weight - nr.weight)/nr.weight

        return R, w

    def _calc_gg_central_satellite(self, data1, data2, one_halo=True, two_halo=True, nbins=20, min_sep=44, max_sep=6e3, slop=0.1, randoms=True, weights=None, return_all=False, verbose=True):

        """This is bloody horrible, I know. But it's the only way I can see to separate out the one and two halo contributions here."""
        
        mask1 = (data1['central']==1)
        mask2 = (data2['central']==0)

        group_ids = np.unique(data1['groupId'])

        if verbose:
            print "Will process:"
            print "1 halo term:", 'yes'*int(one_halo), 'no'*int(np.invert(one_halo))
            print "2 halo term:", 'yes'*int(two_halo), 'no'*int(np.invert(two_halo))
            print ''
            print 'Computing correlations...'

        w1h = []
        w2h = []

        for ig1 in group_ids:
            # This is the central. Should be one galaxy only.
            maski = mask1 & (data1['groupId']==ig1)
            cat_i = treecorr.Catalog(w=weights, x=data1['x'][maski], y=data1['y'][maski], z=data1['z'][maski])
            for ig2 in group_ids:
                if verbose:
                    print ig1,ig2

                if (ig1==ig2) and (not one_halo):
                    continue
                elif (ig1!=ig2) and (not two_halo):
                    continue

                # This is the satellite population for this group.
                maskj = mask2 & (data2['groupId']==ig2)
                cat_j = treecorr.Catalog(w=weights, x=data2['x'][maskj], y=data2['y'][maskj], z=data2['z'][maskj])

                # Build a catalogue of random points drawn from the same volume
                rx_j = np.random.random(size=data2['x'][maskj].size) * (data2['x'][maskj].max()-data1['x'][maskj].min()) + data2['x'][maskj].mean()
                ry_j = np.random.random(size=data2['x'][maskj].size) * (data2['y'][maskj].max()-data1['y'][maskj].min()) + data2['y'][maskj].mean()
                rz_j = np.random.random(size=data2['x'][maskj].size) * (data2['z'][maskj].max()-data1['z'][maskj].min()) + data2['z'][maskj].mean()
                rancat_j  = treecorr.Catalog(x=rx_j, y=ry_j, z=rz_j)

                nn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
                nr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)

                nn.process(cat_i,cat_j) #, metric='Periodic')
                nr.process(rancat_j,cat_j) #, metric='Periodic')

                R = np.exp(nn.meanlogr)
                w = (nn.weight - nr.weight)/nr.weight

                if ig1==ig2:
                    w1h.append([R, w])
                else:
                    w2h.append([R, w])

        return R, w1h, w2h


    def calc_pos_pos(self, i, j, mask1=None, mask2=None, nbins=20, min_sep=44, max_sep=6e3, slop=0.1, randoms=True, weights=None, return_all=False):
        data1 =  self.get(table='subfind_halos' , fields='groupId, x, y, z, mass', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, i))
        data2 =  self.get(table='subfind_halos' , fields='groupId, x, y, z, mass', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, j))
        mask1 = parse_mask(mask1, data1)
        mask2 = parse_mask(mask2, data2)

        # Build a catalogue of random points drawn from the same volume
        rx = np.random.random(size=data1['x'].size) * (data1['x'].max()-data1['x'].min()) + data1['x'].mean()
        ry = np.random.random(size=data1['x'].size) * (data1['y'].max()-data1['y'].min()) + data1['y'].mean()
        rz = np.random.random(size=data1['x'].size) * (data1['z'].max()-data1['z'].min()) + data1['z'].mean()

        # Create the catalogues
        cat_i = treecorr.Catalog(w=weights, x=data1['x'][mask1], y=data1['y'][mask1], z=data1['z'][mask1])
        cat_j = treecorr.Catalog(w=weights, x=data2['x'][mask2], y=data2['y'][mask2], z=data2['z'][mask2])
        rancat_i  = treecorr.Catalog(x=rx, y=ry, z=rz)
        rancat_j  = treecorr.Catalog(x=rx, y=ry, z=rz)

        nn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
        rn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
        nr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
        rr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
        nn.process(cat_i,cat_j) #, metric='Periodic')
        rn.process(rancat_i,cat_j) #, metric='Periodic')
        nr.process(cat_i,rancat_j) #, metric='Periodic')
        rr.process(rancat_i,rancat_j) #, metric='Periodic')

        R = np.exp(nn.meanlogr)
        if randoms:
            w, werr = nn.calculateXi(rr,dr=nr,rd=rn)
        else:
            w, werr = nn.calculateXi(rr,dr=None,rd=None)
        werr = np.sqrt(werr)

        if return_all:
            return R, w, werr, (nn, rn, nr, rr)

        return R, w, werr

    def correlate(self,group1=0,group2=0):
        print 'Will correlate galaxies in groups %d and %d'%(group1,group2)

        data1 =  self.get(table='subfind_halos' , fields='groupId, x, y, z', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, group1))
        data2 =  self.get(table='subfind_halos' , fields='groupId, x, y, z', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, group2))

        #Setup the correlation
        print 'Setting up correlation'
        corr = treecorr.NNCorrelation(nbins=15, min_sep=30, max_sep=4e3)
        cat1 = treecorr.Catalog(x=data1['x'],y=data1['y'],z=data1['z'])
        cat2 = treecorr.Catalog(x=data2['x'],y=data2['y'],z=data2['z'])
        print 'Calculating...'
        corr.process(cat1,cat2)

        print 'Random-Random'
        rr = treecorr.NNCorrelation(nbins=15, min_sep=30, max_sep=4e3)
        rx,ry,rz = np.random.choice(data1['x'], size=5000),np.random.choice(data1['y'], size=5000),np.random.choice(data1['z'], size=5000)
        rcat = treecorr.Catalog(x=rx, y=ry, z=rz)
        rr.process(rcat)

        print 'Data-Random'
        dr = treecorr.NNCorrelation(nbins=15, min_sep=30, max_sep=4e3)
        dr.process(rcat, cat2)

        print 'Random-Data'
        rd = treecorr.NNCorrelation(nbins=15, min_sep=30, max_sep=4e3)
        rd.process(cat1, rcat)

        xi,varxi = corr.calculateXi(rr,dr,rd)

        return np.exp(corr.logr), xi, varxi

    def plot_mass_function(self, groups=[], type='halo'):
        if len(groups)<1:
            groups = self.groups
            self.ngroup = len(self.groups)

        colours = ['purple', 'forestgreen', 'royalblue', 'pink', 'plum', 'darkred']
        ngroup = len(groups)

        for i, group in enumerate(groups):
            print 'Group %d (%d/%d)'%(group, i+1, ngroup)

            if type=='halo':
                name = 'Halo'
                col='mass'
                halo_info = self.get(table='subfind_halos' , fields='groupId, mass', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, group))
                mass_bins = np.linspace(7,13,60)
            elif type=='stellar':
                name = 'Stellar'
                col='m_star'
                halo_info = self.get(table='subfind_halos' , fields='groupId, m_star', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, group))
                mass_bins = np.linspace(-3.5,2,60)
            
            h, bins = np.histogram(np.log10(halo_info[col])[np.isfinite(np.log10(halo_info[col]))], bins=mass_bins, normed=1)


            m = (bins[1:]+bins[:-1])/2
            plt.plot(m, h, color=colours[i], lw=2, label='Group %d'%group)
        
        plt.legend(loc='upper right')
        plt.yscale('log')
        plt.xlabel('%s Mass $\log(M/10^{10} M_\odot h^{-1})$'%name, fontsize=18)
        plt.title("Snapshot %d $z=0.00$"%self.snapshot, fontsize=16)

        print 'Done'


def parse_mask(mask,array):
    if mask is None:
        return np.ones(array.size).astype(bool)
    else:
        return mask


def build_ep_corr(group1, group2, nbins=20, verbose=False):
    rbins = np.linspace(44,6e3, nbins)
    x = (rbins[1:]+rbins[:-1])/2

    y = []

    for lower,upper in zip(rbins[:-1], rbins[1:]):
        vec = []
        print lower, upper
        for i, row1 in enumerate(group1.info):
            for j, row2 in enumerate(group2.info):
                if (row1['subfindId']==row2['subfindId']):
                    if verbose: 
                        print 'Skipping object with subfind ID %d'%row1['subfindId']
                    continue
                R = (row1['x']-row2['x'])**2
                R += (row1['y']-row2['y'])**2
                R += (row1['z']-row2['z'])**2
                R = np.sqrt(R)

                if (R<upper) and (R>lower):
                    a = group1.a[i]
                    r = np.array([(row1['x']-row2['x']), (row1['y']-row2['y']), (row1['z']-row2['z'])])
                    vec.append(np.dot(a,r)*2)
                if verbose:
                    print i, j

        y.append([np.mean(vec), np.std(vec), len(vec)])

    return x, y

def nfw(r,c,a):
    return a/(r/c)/(1+r/c)/(1+r/c)

# Wrapper class that inherits the basic SQL query functions from mbdb
class halos(mbdb):
    def __init__(self, snapshot=85, group=0):
        super(halos, self).__init__()
        self.snapshot=snapshot
        self.group=group

        print "Will load info for halos in group %d snapshot %d"%(group, snapshot)
        return None

    def compile_data(self, nmin=1000):
        data = []
        ngroups = 0
        for j in np.arange(0,10000,1):
            # First submit a query to find the particle membership of this group
            group_info = self.get(table='subfind_groups' , fields='subfindId, nhalo', cond='subfind_groups.snapnum = %d AND subfind_groups.groupId = %d'%(self.snapshot, j))
            nhalo = group_info['nhalo']

            # Check it meets the minimum number of particles
            if nhalo<nmin:
                continue
            # If it does, then get the position data
            particle_info = self.get(table='subfind_halos' , fields='subfindId, groupId, central, mass, len, x, y, z', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, j))

            data.append(particle_info)
            ngroups+=1

        print 'Done'
        print '%d groups'%ngroups
        return np.concatenate(data)

    def get_data(self, info=True, star_shapes=True, dm_shapes=True, star_properties=True):

        if info:
            print 'Loading basic halo information'
            self.info = self.get(table='subfind_halos' , fields='subfindId, central, mass, len, x, y, z, nhalo', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, self.group))
        if star_shapes:
            print 'Loading shapes of stellar halo components'
            self.star_shapes = self.cross_match(self.info, 'subfind_shapes_star', 'subfindId, q2d, a3d_x, a3d_y, a3d_z, b3d_x, b3d_y, b3d_z, c3d_x, c3d_y, c3d_z, a2d_x, a2d_y, b2d_x, b2d_y', 'subfindId', 'subfindId')
            self.a = np.array([self.star_shapes['a3d_x'], self.star_shapes['a3d_y'], self.star_shapes['a3d_z']]).T
            self.b = np.array([self.star_shapes['b3d_x'], self.star_shapes['b3d_y'], self.star_shapes['b3d_z']]).T
            self.c = np.array([self.star_shapes['c3d_x'], self.star_shapes['c3d_y'], self.star_shapes['c3d_z']]).T

        if dm_shapes:
            print 'Loading shapes of dark matter halo components'
            self.dm_shapes= self.cross_match(self.info, 'subfind_shapes_dm', 'subfindId, q2d, a3d_x, a3d_y, a3d_z, b3d_x, b3d_y, b3d_z, c3d_x, c3d_y, c3d_z, a2d_x, a2d_y, b2d_x, b2d_y', 'subfindId', 'subfindId')
        if star_properties:
            print 'Loading extra information on stellar components'
            self.stars = self.cross_match(self.info, 'subfind_star_prop', 'subfindId, btr', 'subfindId')

        print 'Done'
        return None

