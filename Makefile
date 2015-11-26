CGI_DIR=/var/www/bgi-bin/cris4p


all:
	@if [ ! -d $(CGI_DIR) ]; then mkdir $(CGI_DIR);  fi 
	@cp -r * $(CGI_DIR)/
	@echo "crispr4p has been deployed in $(CGI_DIR)"

clean:
	rm -r $(CGI_DIR)/*


