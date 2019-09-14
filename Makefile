all:	growth leveque

growth:
		$(MAKE) -f growth.mk

leveque:
		$(MAKE) -f leveque.mk

clean:
		 rm -rf build/* bin/*
