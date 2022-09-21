all: librte_ecckd.so librte_ecckd.a

librte_ecckd.so: gas_optics_ecckd.o
	$(FC) $(FCFLAGS) -shared -o $@ $^

librte_ecckd.a: gas_optics_ecckd.o
	ar -rvs $@ $^

gas_optics_ecckd.o: src/gas_optics_ecckd.f90
	$(FC) $(FCFLAGS) -fPIC -o $@ -c $<

test: ecckd_rfmip_lw all

ecckd_rfmip_lw: ecckd_rfmip_lw.o \
                mo_load_coefficients.o \
                mo_rfmip_io.o \
                mo_simple_netcdf.o
	$(FC) $(FCFLAGS) -o $@ $(LDFLAGS) $^ -lrte_ecckd $(LIBDIR)/librrtmgp.a $(LIBDIR)/librte.a -lnetcdff -lnetcdf

ecckd_rfmip_lw.o: example/rfmip-rad-irf/ecckd_rfmip_lw.F90 \
                  mo_load_coefficients.o \
                  mo_rfmip_io.o \
                  mo_simple_netcdf.o
	$(FC) $(FCFLAGS) -fPIC -o $@ -c $<

mo_load_coefficients.o: example/rfmip-rad-irf/mo_load_coefficients.F90 \
                        mo_simple_netcdf.o
	$(FC) $(FCFLAGS) -fPIC -o $@ -c $<

mo_rfmip_io.o: example/rfmip-rad-irf/mo_rfmip_io.F90 \
               mo_simple_netcdf.o
	$(FC) $(FCFLAGS) -fPIC -o $@ -c $<

mo_simple_netcdf.o: example/rfmip-rad-irf/mo_simple_netcdf.F90
	$(FC) $(FCFLAGS) -fPIC -o $@ -c $<

clean:
	rm -f gas_optics_ecckd.o librte_ecckd.a
	rm -f ecckd_rfmip_lw ecckd_rfmip_lw.o mo_load_coefficients.o mo_rfmip_io.o mo_simple_netcdf.o
