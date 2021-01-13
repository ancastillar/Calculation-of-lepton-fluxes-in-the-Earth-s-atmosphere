class AIRSAtmosphere(EarthsAtmosphere):

    """Interpolation class for tabulated atmospheres.

    This class is intended to read preprocessed AIRS Satellite data.

    Args:
      location (str): see :func:`init_parameters`
      season (str,optional): see :func:`init_parameters`
    """

    def __init__(self, location, season, extrapolate=True, *args, **kwargs):
        if location != 'SouthPole':
            raise Exception(self.__class__.__name__ +
                "(): Only South Pole location supported. " + location)

        self.extrapolate = extrapolate

        self.month2doy = {'January':1,
                          'February':32,
                          'March':60,
                          'April':91,
                          'May':121,
                          'June':152,
                          'July':182,
                          'August':213,
                          'September':244,
                          'October':274,
                          'November':305,
                          'December':335}
        
        self.hack_list = args[0]
        self.season = season
        
        self.day = self.hack_list['day']
        self.year = self.hack_list['year']
        
        self.init_parameters(location, **kwargs)
        EarthsAtmosphere.__init__(self)

    def init_parameters(self, location, **kwargs):
        """Loads tables and prepares interpolation.

        Args:
          location (str): supported is only "SouthPole"
          doy (int): Day Of Year
        """
        from matplotlib.dates import strpdate2num, num2date
        from os import path

        # YOUR DIRECTORY with subfolder tables/
        data_path = '/data/user/dsoldin/MCEq/MCEq/'
        if 'table_path' in kwargs:
            data_path = kwargs['table_path']

        files = [
            ('dens', 'airs_amsu_dens_180_daily.txt'),
            ('temp', 'airs_amsu_temp_180_daily.txt'),
            ('alti', 'airs_amsu_alti_180_daily.txt')]

        data_collection = {}

        # limit SouthPole pressure to <= 600
        min_press_idx = 4

        IC79_idx_1 = None
        IC79_idx_2 = None

        for d_key, fname in files:
            fname = data_path + 'tables/' + fname
            tab = np.loadtxt(fname,
                             converters={0: strpdate2num('%Y/%m/%d')},
                             usecols=[0] + range(2, 27))
            with open(fname, 'r') as f:
                comline = f.readline()
            # print comline
            p_levels = [float(s.strip()) for s in
                        comline.split(' ')[3:] if s != ''][min_press_idx:]
            dates = num2date(tab[:, 0])
            for di, date in enumerate(dates):
                if date.month == 6 and date.day == 1:
                    if date.year == 2010:
                        IC79_idx_1 = di
                    elif date.year == 2011:
                        IC79_idx_2 = di
            surf_val = tab[:, 1]
            cols = tab[:, min_press_idx + 2:]
            data_collection[d_key] = (dates, surf_val, cols)

        self.interp_tab_d = {}
        self.interp_tab_t = {}
        self.dates = {}
        dates = data_collection['alti'][0]

        msis = MSIS00Atmosphere(location,'January')
        for didx, date in enumerate(dates):
            h_vec = np.array(data_collection['alti'][2][didx,:]*1e2)
            d_vec = np.array(data_collection['dens'][2][didx,:])
            t_vec = np.array(data_collection['temp'][2][didx,:])

            if self.extrapolate:
                #Extrapolate using msis
                h_extra = np.linspace(h_vec[-1],config['h_atm']*1e2,250)
                msis._msis.set_doy(self._get_y_doy(date)[1]-1)
                msis_extra_d = np.array([msis.get_density(h) for h in h_extra])
                msis_extra_t = np.array([msis.get_temperature(h) for h in h_extra])

                # Interpolate last few altitude bins
                ninterp = 5

                for ni in range(ninterp):
                    cl = (1 - np.exp(-ninterp+ni + 1))
                    ch = (1 - np.exp(-ni))
                    norm = 1./(cl + ch)
                    d_vec[-ni-1] = (d_vec[-ni-1]*cl*norm +
                                    msis.get_density(h_vec[-ni-1])*ch*norm)
                    t_vec[-ni-1] = (t_vec[-ni-1]*cl*norm +
                                    msis.get_temperature(h_vec[-ni-1])*ch*norm)

                # Merge the two datasets
                h_vec = np.hstack([h_vec[:-1], h_extra])
                d_vec = np.hstack([d_vec[:-1], msis_extra_d])
                t_vec = np.hstack([t_vec[:-1], msis_extra_t])

            self.interp_tab_d[self._get_y_doy(date)] = (h_vec, d_vec)
            self.interp_tab_t[self._get_y_doy(date)] = (h_vec, t_vec)

            self.dates[self._get_y_doy(date)] = date

        self.IC79_start = self._get_y_doy(dates[IC79_idx_1])
        self.IC79_end = self._get_y_doy(dates[IC79_idx_2])
        self.IC79_days = (dates[IC79_idx_2] - dates[IC79_idx_1]).days
        self.location = location
        if self.season is None:
            self.set_IC79_day(0)
        else:
            self.set_season(self.season)
        # Clear cached value to force spline recalculation
        self.theta_deg = None

    def set_date(self, year, doy):
        self.h, self.dens = self.interp_tab_d[(year, doy)]
        _, self.temp = self.interp_tab_t[(year, doy)]
        self.date = self.dates[(year, doy)]
        # Compatibility with caching
        self.season = self.date

    def _set_doy(self, doy, year=2010):
        self.h, self.dens = self.interp_tab_d[(year, doy)]
        _, self.temp = self.interp_tab_t[(year, doy)]
        self.date = self.dates[(year, doy)]
        
        print self.date

    def set_season(self, month):
        self.season = month
        if self.day==0:
	        self._set_doy(self.month2doy[month])
        else:
        	self._set_doy(self.day, self.year)
        self.season = month

    def set_IC79_day(self, IC79_day):
        import datetime
        if IC79_day > self.IC79_days:
            raise Exception(self.__class__.__name__ +
                            "::set_IC79_day(): IC79_day above range.")
        target_day = self._get_y_doy(self.dates[self.IC79_start] +
                                     datetime.timedelta(days=IC79_day))
        if dbg:
            print 'setting IC79_day', IC79_day
        self.h, self.dens = self.interp_tab_d[target_day]
        _, self.temp = self.interp_tab_t[target_day]
        self.date = self.dates[target_day]
        # Compatibility with caching
        self.season = self.date

    def _get_y_doy(self, date):
        return date.timetuple().tm_year, date.timetuple().tm_yday

    def get_density(self, h_cm):
        """ Returns the density of air in g/cm**3.

        Interpolates table at requested value for previously set
        year and day of year (doy).

        Args:
          h_cm (float): height in cm

        Returns:
          float: density :math:`\\rho(h_{cm})` in g/cm**3
        """
        ret = np.exp(np.interp(h_cm, self.h, np.log(self.dens)))
        try:
            ret[h_cm > self.h[-1]] = np.nan
        except TypeError:
            if h_cm > self.h[-1]: return np.nan
        return ret

    def get_temperature(self, h_cm):
        """ Returns the temperature in K.

        Interpolates table at requested value for previously set
        year and day of year (doy).

        Args:
          h_cm (float): height in cm

        Returns:
          float: temperature :math:`T(h_{cm})` in K
        """
        ret = np.exp(np.interp(h_cm, self.h, np.log(self.temp)))
        try:
            ret[h_cm > self.h[-1]] = np.nan
        except TypeError:
            if h_cm > self.h[-1]: return np.nan
        return ret