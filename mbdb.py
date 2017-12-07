#import MySQLdb as mdb
import pymysql as mdb
import numpy as np
from numpy.core.records import fromarrays
import pylab as plt
import treecorr
#import tools.plots as pl
plt.switch_backend('pdf')
import halotools as ht
import halotools.mock_observables as pretending
import mbii.lego_tools as util
import mbii.basic_simulation_info as info
import tools.diagnostics as di
import fitsio as fi

#plt.style.use('y1a1')


class mbdb(object):

    def __init__(self, sqlserver='localhost', user='flanusse', password='mysqlpass',
                 dbname='mb2_hydro', unix_socket='/home/rmandelb.proj/flanusse/mysql/mysql.sock', fatal_errors=True):
        try:
            self.db = mdb.connect(sqlserver, user, password, dbname, unix_socket=unix_socket)
        except:
            if fatal_errors:
                raise 
            else:
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
            sql+="'%d',"%int(row)
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

        # Finally match the results
        # Without this the results of the second query are misaligned  
        sm, rm = di.match_results(source, results, name1='subfindId', name2='subfindId')
        return sm, rm

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

class particle_data(mbdb):
    def __init__(self, filename):
        super(particle_data, self).__init__()
        print 'Obtaining particle catalogue from'
        print filename
        self.info = fi.FITS(filename)['particle_data'].read()

class groups(mbdb):
    def __init__(self, snapshot=85, info=None):
        super(groups, self).__init__()

        print "Will load group data for snapshot %d"%(snapshot)
        self.snapshot=snapshot
        self.Lbox = 100 * 1e3

        if info is not None:
        	self.info = info

        #info = self.get(table='subfind_halos' , fields='groupId', cond='subfind_halos.snapnum = %d'%(self.snapshot))
        #self.groups = np.unique(info['groupId'])
        return None

    def _get_corrs_nosep(self, data, min_sep=44, max_sep=1e6, binning='log', nbins=20, ctype=('s','s'), estimator='Landy-Szalay', verbosity=1, randoms=None, method='halotools'):

    	if verbosity>0:
    		print 'Will construct %s - %s correlation functions'%ctype
    		print 'Using %s estimator'%estimator

    	# Decide on an appropriate binning scheme
    	if (binning.lower()=='log'):
    		rbins = np.logspace(np.log10(min_sep), np.log10(max_sep), nbins )
    	elif (binning.lower()=='linear'):
    		rbins = np.linspace(min_sep, max_sep, nbins )

    	if verbosity>1:
    		print 'Will use %s binning:'%binning, rbins

    	# Parse the mask
    	mask1 = util.choose_cs_mask(data,ctype[0])
    	mask2 = util.choose_cs_mask(data,ctype[1])

    	pos1 = pretending.return_xyz_formatted_array(data['x'], data['y'], data['z'], mask = mask1)
    	pos2 = pretending.return_xyz_formatted_array(data['x'], data['y'], data['z'], mask = mask2)

    	# And do the randoms
    	if randoms is None:
            r1 = util.construct_random_cat(data, mask=mask1)
            r2 = util.construct_random_cat(data, mask=mask2)
        else:
            if verbosity>0:
                print 'Using random points provided for normalisation.'
            r1 = randoms

    	R = np.sqrt(np.array(rbins)[1:]*np.array(rbins)[:-1]) 

        print 'Using %s to calculate two-point correlations'%method
    	if method=='halotools':
            return R, pretending.tpcf(pos1, rbins, sample2=pos2, randoms=r1, period=info.Lbox, estimator=estimator )

        elif method=='treecorr':
            print 'Constructing catalogues...'
            cat_i = treecorr.Catalog(x=data['x'][mask1], y=data['y'][mask1], z=data['z'][mask1])
            cat_j = treecorr.Catalog(x=data['x'][mask2], y=data['y'][mask2], z=data['z'][mask2])
            rx_1 = (np.random.random(size=data['x'][mask1].size) - 0.5) * (data['x'][mask1].max()-data['x'][mask1].min()) + data['x'][mask1].mean()
            ry_1 = (np.random.random(size=data['x'][mask1].size) - 0.5) * (data['y'][mask1].max()-data['y'][mask1].min()) + data['y'][mask1].mean()
            rz_1 = (np.random.random(size=data['x'][mask1].size) - 0.5) * (data['z'][mask1].max()-data['z'][mask1].min()) + data['z'][mask1].mean()
            rancat_1  = treecorr.Catalog(x=rx_1, y=ry_1, z=rz_1)

            print 'Correlating...'
            nn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=0.1)
            nr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=0.1)
            rn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=0.1)
            rr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=0.1)

            nn.process(cat_i,cat_j) 
            nr.process(rancat_1,cat_i)
            rn.process(cat_j,rancat_1)
            rr.process(rancat_1,rancat_1) 
            R = np.exp(nn.meanlogr)
            w = (nn.weight - nr.weight -  rn.weight + rr.weight) / rr.weight
            return R, w



    def _get_corrs(self, data, min_sep=44, max_sep=1e6, binning='log', nbins=20, ctype=('s','s'), estimator='Landy-Szalay', verbosity=1, fran=(1,1), randoms=None):

    	if verbosity>0:
    		print 'Will construct %s - %s correlation functions'%ctype
    		print 'Using %s estimator'%estimator

    	# Decide on an appropriate binning scheme
    	if (binning.lower()=='log'):
    		rbins = np.logspace(np.log10(min_sep), np.log10(max_sep), nbins )
    	elif (binning.lower()=='linear'):
    		rbins = np.linspace(min_sep, max_sep, nbins )

    	if verbosity>1:
    		print 'Will use %s binning:'%binning, rbins

    	# Parse the mask
    	mask1 = util.choose_cs_mask(data,ctype[0])
    	mask2 = util.choose_cs_mask(data,ctype[1])

    	pos1 = pretending.return_xyz_formatted_array(data['x'], data['y'], data['z'], mask = mask1)
    	pos2 = pretending.return_xyz_formatted_array(data['x'], data['y'], data['z'], mask = mask2)

    	# And do the randoms
    	if randoms is None:
            r1 = util.construct_random_cat(data, mask=mask1, f=fran[0])
            r2 = util.construct_random_cat(data, mask=mask2, f=fran[1])
        else:
            if verbosity>0:
                print 'Using random points provided for normalisation.'
            r1 = randoms

    	R = np.sqrt(np.array(rbins)[1:]*np.array(rbins)[:-1]) 

        #print 'Computing jackkife errorbars'
        #E11 = pretending.tpcf_jackknife(pos1, r1, rbins, sample2=pos1, period=self.Lbox, estimator=estimator )
        #E22 = pretending.tpcf_jackknife(pos2, r1, rbins, sample2=pos2, period=self.Lbox, estimator=estimator )
        #E12 = pretending.tpcf_jackknife(pos1, r1, rbins, sample2=pos2, period=self.Lbox, estimator=estimator )

    	return R, pretending.tpcf_one_two_halo_decomp(pos1, data['groupId'][mask1], rbins, sample2=pos2,  sample2_host_halo_id=data['groupId'][mask2], randoms=r1, randoms2=None, factor=fran, period=self.Lbox, estimator=estimator )

    	#return rbins, xi_1h_11, xi_2h_11, xi_1h_12, xi_2h_12, xi_1h_22, xi_2h_22



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

    def correlate_all(self, data1, data2, pair=('c','s'), halos=2, nbins=20, slop=0.1, min_sep=44, max_sep=6e3):

        #mask1 = (data1['central']==int(pair[0]=='c'))
        #mask2 = (data2['central']==int(pair[1]=='c'))

        rx_j = (np.random.random(size=data2['x'].size) - 0.5)* (data2['x'].max()-data2['x'].min()) + data2['x'].mean()
        ry_j = (np.random.random(size=data2['x'].size) - 0.5)* (data2['y'].max()-data2['y'].min()) + data2['y'].mean()
        rz_j = (np.random.random(size=data2['x'].size) - 0.5)* (data2['z'].max()-data2['z'].min()) + data2['z'].mean()
        rancat_i  = treecorr.Catalog(x=rx_j, y=ry_j, z=rz_j)

        rx_i = (np.random.random(size=data1['x'].size) - 0.5) * (data1['x'].max()-data1['x'].min()) + data1['x'].mean()
        ry_i = (np.random.random(size=data1['x'].size) - 0.5) * (data1['y'].max()-data1['y'].min()) + data1['y'].mean()
        rz_i = (np.random.random(size=data1['x'].size) - 0.5) * (data1['z'].max()-data1['z'].min()) + data1['z'].mean()
        rancat_i  = treecorr.Catalog(x=rx_i, y=ry_i, z=rz_i)

        nn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
        nr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
        rn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
        rr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)

        nn.process(cat_i,cat_j) #, metric='Periodic')
        nr.process(rancat_j,cat_i) #, metric='Periodic')
        rn.process(cat_j,rancat_i) #, metric='Periodic')
        rr.process(rancat_i,rancat_j) #, metric='Periodic')

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

    def _calc_2h_cs(self, data1, data2, mask1, mask2, save=False, verbose=False, weights=None, nbins=20, min_sep=44, max_sep=6e3, slop=0.1):
        """Given two numpy arrays of positions, compute the two halo central-satellite realspace correlation."""

        w2h_cs = []
        group_ids = np.unique(data1['groupId'])
        N = len(group_ids)

        for ig1 in group_ids:
            if verbose: 
                print '%d/%d'%(ig1+1, N)
            maski = mask1 & (data1['groupId']==ig1)
            cat_i = treecorr.Catalog(w=weights, x=data1['x'][maski], y=data1['y'][maski], z=data1['z'][maski])

            maskj = mask2 & (data2['groupId']!=ig1)
            cat_j = treecorr.Catalog(w=weights, x=data2['x'][maskj], y=data2['y'][maskj], z=data2['z'][maskj])

            rx_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['x'][maskj].max()-data2['x'][maskj].min()) + data2['x'][maskj].mean()
            ry_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['y'][maskj].max()-data2['y'][maskj].min()) + data2['y'][maskj].mean()
            rz_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['z'][maskj].max()-data2['z'][maskj].min()) + data2['z'][maskj].mean()
            rancat_j  = treecorr.Catalog(x=rx_j, y=ry_j, z=rz_j)
                
            f=10000
            rx_i = (np.random.random(size=data1['x'][maski].size * f) - 0.5) * (data2['x'][maskj].max()-data2['x'][maskj].min()) + data2['x'][maskj].mean()
            ry_i = (np.random.random(size=data1['x'][maski].size * f) - 0.5) * (data2['y'][maskj].max()-data2['y'][maskj].min()) + data2['y'][maskj].mean()
            rz_i = (np.random.random(size=data1['x'][maski].size * f) - 0.5) * (data2['z'][maskj].max()-data2['z'][maskj].min()) + data2['z'][maskj].mean()
            rancat_i  = treecorr.Catalog(x=rx_i, y=ry_i, z=rz_i)

            nn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            nr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            rn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            rr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)

            nn.process(cat_i,cat_j) #, metric='Periodic')
            nr.process(rancat_j,cat_i) #, metric='Periodic')
            rn.process(cat_j,rancat_i) #, metric='Periodic')
            rr.process(rancat_i,rancat_j) #, metric='Periodic')

            R_2h_cs = np.exp(nn.meanlogr)
            coeff = 1./f
            w = (nn.weight - nr.weight - coeff * rn.weight + coeff*rr.weight)/(coeff * rr.weight)
            w2h_cs.append(w)

        if save:
            print 'Storing...'
            np.savetxt('R_2h_cs.txt', R_2h_cs)
            np.savetxt('w2h_cs.txt', w2h_cs)

        return R_2h_cs, w2h_cs

    def _calc_2h_cc(self, data1, data2, mask1, mask2, save=False, verbose=False, weights=None, nbins=20, min_sep=44, max_sep=6e3, slop=0.1):

        w2h_cc=[]
        group_ids = np.unique(data1['groupId'])
        N = len(group_ids)

        for ig1 in group_ids:
            if verbose: 
                print '%d/%d'%(ig1+1, N)
            maski = mask1 & (data1['groupId']==ig1)
            cat_i = treecorr.Catalog(w=weights, x=data1['x'][maski], y=data1['y'][maski], z=data1['z'][maski])

            # Select all of the centrals that are not part of the same halo
            maskj = mask1 & (data1['groupId']!=ig1)
            cat_j = treecorr.Catalog(w=weights, x=data2['x'][maskj], y=data2['y'][maskj], z=data2['z'][maskj])

            rx_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['x'][maskj].max()-data2['x'][maskj].min()) + data2['x'][maskj].mean()
            ry_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['y'][maskj].max()-data2['y'][maskj].min()) + data2['y'][maskj].mean()
            rz_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['z'][maskj].max()-data2['z'][maskj].min()) + data2['z'][maskj].mean()
            rancat_j  = treecorr.Catalog(x=rx_j, y=ry_j, z=rz_j)
                
            f=10000
            rx_i = (np.random.random(size=data1['x'][maski].size * f) -0.5) * (data2['x'][maskj].max()-data2['x'][maskj].min()) + data2['x'][maskj].mean()
            ry_i = (np.random.random(size=data1['x'][maski].size * f) -0.5) * (data2['y'][maskj].max()-data2['y'][maskj].min()) + data2['y'][maskj].mean()
            rz_i = (np.random.random(size=data1['x'][maski].size * f) -0.5) * (data2['z'][maskj].max()-data2['z'][maskj].min()) + data2['z'][maskj].mean()
            rancat_i  = treecorr.Catalog(x=rx_i, y=ry_i, z=rz_i)

            nn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            nr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            rn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            rr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)

            nn.process(cat_i,cat_j) #, metric='Periodic')
            nr.process(rancat_j,cat_i) #, metric='Periodic')
            rn.process(cat_j,rancat_i) #, metric='Periodic')
            rr.process(rancat_i,rancat_j) #, metric='Periodic')

            R_2h_cc = np.exp(nn.meanlogr)
            coeff = 1./f
            w = (nn.weight - nr.weight - coeff * rn.weight + coeff*rr.weight)/(coeff * rr.weight)
            w2h_cc.append(w)

        if save:
            print 'Storing...'
            np.savetxt('R_2h_cc.txt',R_2h_cc)
            np.savetxt('w2h_cc.txt',w2h_cc)

        return R_2h_cc, w2h_cc

    def _calc_2h_ss(self, data1, data2, mask1, mask2, save=False, verbose=False, weights=None, nbins=20, min_sep=44, max_sep=6e3, slop=0.1):

        w2h_ss=[]
        group_ids = np.unique(data1['groupId'])
        N = len(group_ids)

        for j, ig1 in enumerate(group_ids):
            if verbose: 
                print '%d %d/%d'%(ig1, j+1, N)
            maski = mask2 & (data1['groupId']==ig1)
            cat_i = treecorr.Catalog(w=weights, x=data1['x'][maski], y=data1['y'][maski], z=data1['z'][maski])

            maskj = mask2 & (data2['groupId']!=ig1)
            cat_j = treecorr.Catalog(w=weights, x=data2['x'][maskj], y=data2['y'][maskj], z=data2['z'][maskj])

            rx_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['x'][maskj].max()-data2['x'][maskj].min()) + data2['x'][maskj].mean()
            ry_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['y'][maskj].max()-data2['y'][maskj].min()) + data2['y'][maskj].mean()
            rz_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['z'][maskj].max()-data2['z'][maskj].min()) + data2['z'][maskj].mean()
            rancat_j  = treecorr.Catalog(x=rx_j, y=ry_j, z=rz_j)

            rx_i = (np.random.random(size=data1['x'][maski].size) - 0.5) * (data2['x'][maskj].max()-data2['x'][maskj].min()) + data2['x'][maskj].mean()
            ry_i = (np.random.random(size=data1['x'][maski].size) - 0.5) * (data2['y'][maskj].max()-data2['y'][maskj].min()) + data2['y'][maskj].mean()
            rz_i = (np.random.random(size=data1['x'][maski].size) - 0.5) * (data2['z'][maskj].max()-data2['z'][maskj].min()) + data2['z'][maskj].mean()
            rancat_i  = treecorr.Catalog(x=rx_i, y=ry_i, z=rz_i)

            nn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            nr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            rn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            rr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)

            nn.process(cat_i,cat_j) #, metric='Periodic')
            nr.process(rancat_j,cat_i) #, metric='Periodic')
            rn.process(cat_j,rancat_i) #, metric='Periodic')
            rr.process(rancat_i,rancat_j) #, metric='Periodic')

            R_2h_ss = np.exp(nn.meanlogr)
            w = (nn.weight - nr.weight -  rn.weight + rr.weight) / rr.weight
            w2h_ss.append(w)

        if save:
            print 'Storing...'
            np.savetxt('R_2h_ss.txt', R_2h_ss)
            np.savetxt('w2h_ss.txt', w2h_ss)

        return R_2h_ss, w2h_ss

    def _calc_1h_cs(self, data1, data2, mask1, mask2, save=False, verbose=False, weights=None, nbins=20, min_sep=44, max_sep=6e3, slop=0.1):

        w1h_cs=[]
        group_ids = np.unique(data1['groupId'])
        N = len(group_ids)

        for ig1 in group_ids:
            if verbose: 
                print '%d/%d'%(ig1+1, N)
            maski = mask1 & (data1['groupId']==ig1)
            cat_i = treecorr.Catalog(w=weights, x=data1['x'][maski], y=data1['y'][maski], z=data1['z'][maski])

            maskj = mask2 & (data2['groupId']==ig1)
            cat_j = treecorr.Catalog(w=weights, x=data2['x'][maskj], y=data2['y'][maskj], z=data2['z'][maskj])

            rx_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['x'][maskj].max()-data2['x'][maskj].min()) + data2['x'][maskj].mean()
            ry_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['y'][maskj].max()-data2['y'][maskj].min()) + data2['y'][maskj].mean()
            rz_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['z'][maskj].max()-data2['z'][maskj].min()) + data2['z'][maskj].mean()
            rancat_j  = treecorr.Catalog(x=rx_j, y=ry_j, z=rz_j)
                
            f=10000
            rx_i = (np.random.random(size=data1['x'][maski].size * f) -0.5) * (data2['x'][maskj].max()-data2['x'][maskj].min()) + data2['x'][maskj].mean()
            ry_i = (np.random.random(size=data1['x'][maski].size * f) -0.5) * (data2['y'][maskj].max()-data2['y'][maskj].min()) + data2['y'][maskj].mean()
            rz_i = (np.random.random(size=data1['x'][maski].size * f) -0.5) * (data2['z'][maskj].max()-data2['z'][maskj].min()) + data2['z'][maskj].mean()
            rancat_i  = treecorr.Catalog(x=rx_i, y=ry_i, z=rz_i)

            nn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            nr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            rn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
            rr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)

            nn.process(cat_i,cat_j) #, metric='Periodic')
            nr.process(rancat_j,cat_i) #, metric='Periodic')
            rn.process(cat_j,rancat_i) #, metric='Periodic')
            rr.process(rancat_i,rancat_j) #, metric='Periodic')

            R_1h_cs = np.exp(nn.meanlogr)
            coeff = 1./f
            w = (nn.weight - nr.weight - coeff * rn.weight + coeff*rr.weight)/(coeff * rr.weight)
            w1h_cs.append(w)

        if save:
            print 'Storing...'
            np.savetxt('R_1h_cs.txt',R_1h_cs)
            np.savetxt('w1h_cs.txt',w1h_cs)

        return R_1h_cs, w1h_cs

    def calc_gg_all(self, data1, data2, one_halo=True, corrs=['cc','ss','cs'], two_halo=True, nbins=20, min_sep=44, max_sep=6e3, slop=0.1, randoms=True, weights=None, return_all=False, verbose=True, bootstrap_errors=False, save=False):
        """This is bloody horrible, I know. It's currently the only way I can see to separate out the one and two halo contributions here."""
        
        mask1 = (data1['central']==1)
        mask2 = (data2['central']==0)

        group_ids = np.unique(data1['groupId'])

        if verbose:
            print "Will process:"
            print "1 halo term:", 'yes'*int(one_halo), 'no'*int(np.invert(one_halo))
            print "2 halo term:", 'yes'*int(two_halo), 'no'*int(np.invert(two_halo))
            print ''
            if save:
                print 'Will save text output.'
                print ''
            print 'Computing correlations...'

        R_2h_cs =[]
        w2h_cs =[]
        R_2h_cc =[]
        w2h_cc =[]
        R_2h_ss =[]
        w2h_ss =[]
        R_1h_cs =[]
        w1h_cs =[]

        # Two halo, central - satellite
        if verbose:
            print 'Two halo, central - satellite'
        if (two_halo) and ('cs' in corrs):
            R_2h_cs, w2h_cs = self._calc_2h_cs(data1,data2,mask1,mask2, save=save, verbose=verbose, weights=weights, nbins=nbins, min_sep=min_sep, max_sep=max_sep)

        # Two halo, central - central
        if verbose:
            print 'Two halo, central - central'
        if (two_halo) and ('cc' in corrs):
            R_2h_cc, w2h_cc = self._calc_2h_cc(data1,data2,mask1,mask2, save=save, verbose=verbose, weights=weights, nbins=nbins, min_sep=min_sep, max_sep=max_sep)

        # Two halo, satellite - satellite
        if verbose:
            print 'Two halo, satellite - satellite'
        if (two_halo) and ('ss' in corrs):
            R_2h_ss, w2h_ss = self._calc_2h_ss(data1,data2,mask1,mask2, save=save, verbose=verbose, weights=weights, nbins=nbins, min_sep=min_sep, max_sep=max_sep)

        # One halo, central - satellite
        if verbose:
            print 'One halo, central - satellite'
        if (one_halo) and ('cs' in corrs):
            R_1h_cs, w1h_cs = self._calc_1h_cs(data1,data2,mask1,mask2, save=save, verbose=verbose, weights=weights, nbins=nbins, min_sep=min_sep, max_sep=max_sep)

        return (R_1h_cs, R_2h_cs, R_2h_ss, R_2h_cc), (w1h_cs, w2h_cs, w2h_ss, w2h_cc)

    def _calc_gg_central_satellite_discrete(self, data1, data2, one_halo=True, two_halo=True, nbins=20, min_sep=44, max_sep=6e3, slop=0.1, randoms=True, weights=None, return_all=False, verbose=True, bootstrap_errors=False):
        """This is bloody horrible, I know. It's currently the only way I can see to separate out the one and two halo contributions here."""
        
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

                #import pdb ; pdb.set_trace()

                # Build a catalogue of random points drawn from the same volume
                rx_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['x'][maskj].max()-data2['x'][maskj].min()) + data2['x'][maskj].mean()
                ry_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['y'][maskj].max()-data2['y'][maskj].min()) + data2['y'][maskj].mean()
                rz_j = (np.random.random(size=data2['x'][maskj].size) - 0.5) * (data2['z'][maskj].max()-data2['z'][maskj].min()) + data2['z'][maskj].mean()
                rancat_j  = treecorr.Catalog(x=rx_j, y=ry_j, z=rz_j)

                f=10000
                rx_i = (np.random.random(size=data1['x'][maski].size * f) -0.5) * (data2['x'][maskj].max()-data2['x'][maskj].min()) + data2['x'][maskj].mean()
                ry_i = (np.random.random(size=data1['x'][maski].size * f) -0.5) * (data2['y'][maskj].max()-data2['y'][maskj].min()) + data2['y'][maskj].mean()
                rz_i = (np.random.random(size=data1['x'][maski].size * f) -0.5) * (data2['z'][maskj].max()-data2['z'][maskj].min()) + data2['z'][maskj].mean()
                rancat_i  = treecorr.Catalog(x=rx_i, y=ry_i, z=rz_i)

                nn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
                nr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
                rn = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)
                rr = treecorr.NNCorrelation(nbins=nbins, min_sep=min_sep, max_sep=max_sep, bin_slop=slop)

                nn.process(cat_i,cat_j) #, metric='Periodic')
                nr.process(rancat_j,cat_i) #, metric='Periodic')
                rn.process(cat_j,rancat_i) #, metric='Periodic')
                rr.process(rancat_i,rancat_j) #, metric='Periodic')

                #if bootstrap_errors:
                #    Ew = self.bootstrap_gg_errors(data1, data2, maski, maskj)

                #import pdb ; pdb.set_trace()

                R = np.exp(nn.meanlogr)
                coeff = 1./f
                w = (nn.weight - nr.weight - coeff * rn.weight + coeff*rr.weight)/(coeff * rr.weight)

                if ig1==ig2:
                    w1h.append(w)
                else:
                    w2h.append(w)

        return R, w1h, w2h

    def bootstrap_gg_errors(self, data1, data2, maski, maskj, npatch=20):

        #Decide how many particles to use per patch
        n = data1['x'][maski].size/npatch

        for i in xrange(npatch):
            # Select a random subset of particles
            indices = np.random.choice(data1['x'][maski].size, size=n, replace=False)
            if len(data1['x'][maski])>0:
                cat_i = treecorr.Catalog(w=weights, x=data1['x'][maski], y=data1['y'][maski], z=data1['z'][maski])
            else:
                cat_i = treecorr.Catalog(w=weights, x=data1['x'][maski], y=data1['y'][maski], z=data1['z'][maski])

                cat_j = treecorr.Catalog(w=weights, x=data2['x'][maskj], y=data2['y'][maskj], z=data2['z'][maskj])


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

# Wrapper class that inherits the basic SQL query functions from mbdb
class halos(mbdb):
    def __init__(self, snapshot=85, group=0, fatal_errors=True):
        super(halos, self).__init__()
        self.snapshot=snapshot
        self.group=group

        print "Will load info for halos in group %d snapshot %d"%(group, snapshot)
        return None

    def compile_data(self, nmin=1000, dm_shapes=True, star_shapes=True):
        data = []
        ngroups = 0
        for j in np.arange(0,20000,1):
            # First submit a query to find the particle membership of this group
            group_info = self.get(table='subfind_groups' , fields='subfindId, nhalo', cond='subfind_groups.snapnum = %d AND subfind_groups.groupId = %d'%(self.snapshot, j))
            nhalo = group_info['nhalo']

            # Check it has the minimum number of particles
            if (nhalo<nmin):
                continue
            # If it does, then get the position data
            particle_info = self.get(table='subfind_halos' , fields='subfindId, groupId, central, mass, m_dm, m_star, m_gas, len, x, y, z', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, j))
            if star_shapes:
                names = 'subfindId, q2d, q3d, s3d, a3d_x, a3d_y, a3d_z, b3d_x, b3d_y, b3d_z, c3d_x, c3d_y, c3d_z, a2d_x, a2d_y, b2d_x, b2d_y'
                particle_info, star_info = self.cross_match(particle_info, 'subfind_shapes_star', names, 'subfindId', 'subfindId')
                for col in star_info.dtype.names:
                    if col=='subfindId':
                        continue
                    particle_info = util.add_col(particle_info, '%s_star'%col, star_info[col] )
            if dm_shapes:
                names = 'subfindId, q2d, q3d, s3d, a3d_x, a3d_y, a3d_z, b3d_x, b3d_y, b3d_z, c3d_x, c3d_y, c3d_z, a2d_x, a2d_y, b2d_x, b2d_y'
                particle_info, dm_info = self.cross_match(particle_info, 'subfind_shapes_dm', names, 'subfindId', 'subfindId')
                for col in dm_info.dtype.names:
                    if col=='subfindId':
                        continue
                    particle_info = util.add_col(particle_info, '%s_dm'%col, dm_info[col] )

                # Check Ananth has generated a shape for everything selected 
                # We're using the major axis here. The original cut was on the 'lenbytype' field,
                # which isn't included in the coma sql tables.
                select = (particle_info['a3d_x_dm']!=0) & (particle_info['a3d_y_dm']!=0)
                print 'Discarding %d/%d particles with len<50'%(particle_info['subfindId'][np.invert(select)].size, particle_info['subfindId'].size)
                particle_info = particle_info[select]

            data.append(particle_info)
            ngroups+=1

        print 'Done'
        print '%d groups'%ngroups
        return np.concatenate(data)

    def get_data(self, info=True, star_shapes=True, dm_shapes=True, star_properties=True):

        if info:
            print 'Loading basic halo information'
            self.info = self.get(table='subfind_halos' , fields='subfindId, groupId, central, mass, len, x, y, z', cond='subfind_halos.snapnum = %d AND subfind_halos.groupId = %d'%(self.snapshot, self.group))
        if star_shapes:
            print 'Loading shapes of stellar halo components'
            self.star_shapes = self.cross_match(self.info, 'subfind_shapes_star', 'subfindId, q2d, a3d_x, a3d_y, a3d_z, b3d_x, b3d_y, b3d_z, c3d_x, c3d_y, c3d_z, a2d_x, a2d_y, b2d_x, b2d_y', 'subfindId', 'subfindId')
            #self.a = np.array([self.star_shapes['a3d_x'], self.star_shapes['a3d_y'], self.star_shapes['a3d_z']]).T
            #self.b = np.array([self.star_shapes['b3d_x'], self.star_shapes['b3d_y'], self.star_shapes['b3d_z']]).T
            #self.c = np.array([self.star_shapes['c3d_x'], self.star_shapes['c3d_y'], self.star_shapes['c3d_z']]).T

        if dm_shapes:
            print 'Loading shapes of dark matter halo components'
            self.dm_shapes= self.cross_match(self.info, 'subfind_shapes_dm', 'subfindId, q2d, a3d_x, a3d_y, a3d_z, b3d_x, b3d_y, b3d_z, c3d_x, c3d_y, c3d_z, a2d_x, a2d_y, b2d_x, b2d_y', 'subfindId', 'subfindId')
        if star_properties:
            print 'Loading extra information on stellar components'
            self.stars = self.cross_match(self.info, 'subfind_star_prop', 'subfindId, btr', 'subfindId')

        print 'Done'
        return None

import math
import fitsio as fi

def classify_subhalos(data):
    # Setup the database connection
    sqlserver='localhost'
    user='flanusse'
    password='mysqlpass'
    dbname='mb2_hydro'
    unix_socket='/home/rmandelb.proj/flanusse/mysql/mysql.sock'
    db = mdb.connect(sqlserver, user, password, dbname, unix_socket=unix_socket)
    names = ['subfindId', 'snapnum', 'haloId', 'groupId', 'central', 'mass', 'len', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'vdisp', 'vcirc', 'rcirc', 'm_gas', 'm_dm', 'm_star', 'm_bh'] 

    dt = [('subfindId',int),
    ('snapnum',int),
    ('haloId',int),
    ('groupId',int),
    ('central',int),
    ('mass',float),
    ('len',int),
    ('x',float),
    ('y',float),
    ('z',float),
    ('vx',float),
    ('vy',float),
    ('vz',float),
    ('vdisp',float),
    ('vcirc',float), 
    ('rcirc',float),
    ('m_gas',float),
    ('m_dm',float),
    ('m_star',float),
    ('m_bh',float)]

    out = np.zeros(len(data), dtype=dt)

    for i, subhalo in enumerate(data):
        sql = "SELECT * FROM subfind_halos WHERE len=%d AND groupId=%d;"%(subhalo['len'], math.ceil(subhalo['groupid']))

        #print sql
        cursor = db.cursor()
        cursor.execute(sql)

        try:
            
            results = fromarrays(np.array(cursor.fetchall()).squeeze().T, names=names)
        except:
            print 'Could not match results'
            continue
        print i, results

        if (results.size>1):
            select = np.isclose(results['x'],subhalo['pos'][0]) & np.isclose(results['y'],subhalo['pos'][1]) & np.isclose(results['z'],subhalo['pos'][2])
            results = results[select]
        
        for colname in results.dtype.names:
            out[colname][i] = results[colname]

    outfits = fi.FITS('/home/ssamurof/subhalo_central_flag.fits','rw')
    outfits.write(out)
    outfits.close()











