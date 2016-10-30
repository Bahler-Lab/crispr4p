CGI_DIR=/var/www/bgi-bin/crispr4p


all:
	@if [ ! -d $(CGI_DIR) ]; then mkdir $(CGI_DIR);  fi 
	@cp -r * $(CGI_DIR)/
	@echo "crispr4p has been deployed in $(CGI_DIR)"
	#solucionar lo del css moverlo a ../../bahlerweb.css

regression:
	crispr4p/regression.py


clean:
	rm -r $(CGI_DIR)/*


