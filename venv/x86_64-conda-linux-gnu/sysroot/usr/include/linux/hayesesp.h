#ifndef HAYESESP_H
#define HAYESESP_H

struct hayes_esp_config {
	short flow_on;
	short flow_off;
	short rx_trigger;
	short tx_trigger;
	short pio_threshold;
	unsigned char rx_timeout;
	char dma_channel;
};



#endif /* ESP_H */

