Clair3/models/ont_guppy5/pileup.index: Clair3/run_clair3.sh workflow/scripts/longphas
	cd Clair3 \
		&& mkdir models \
		&& wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz \
		&& tar -zxvf clair3_models.tar.gz -C ./models \
		&& rm clair3_models.tar.gz \
		&& cd ..

Clair3/run_clair3.sh:
	git clone git@github.com:HKU-BAL/Clair3.git \
		&& cd Clair3 \
		&& git checkout v0.1-r7 \
		&& cd ..

workflow/scripts/longphase:
	wget https://github.com/twolinin/longphase/releases/download/v1.4/longphase_linux-x64.tar.xz \
		&& tar -xJf longphase_linux-x64.tar.xz \
		&& mv longphase_linux-x64 workflow/scripts/longphase

.PHONY: clean

clean:
	rm -rf Clair3
