PROGRAM	= vbgmm
CC	= gcc
CFLAGS	= -O3
SRCS	= vbgmm.c util.c writer.c matrix.c lb.c learn.c
OBJS	= $(SRCS:.c=.o)
HEADERS	= $(SRCS:.c=.h)
LDFLAGS	= -lm -lgsl -lgslcblas
VERSION	= 0.1
PKGNAME	= vbgmm-$(VERSION)
DISTDIR	= ../dist
DISTFILES	= $(SRCS) $(HEADERS) Makefile

all: depend $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)
.c.o:
	$(CC) $(CFLAGS) -c $<

depend:
	@$(CC) -MM $(SRCS) > .depend
clean:
	@rm -f .depend $(OBJS)
pkg:
	@[ -d $(PKGNAME) ] || mkdir $(PKGNAME)
	@cp -p $(DISTFILES) $(PKGNAME)
	@tar czvf $(DISTDIR)/$(PKGNAME).tar.gz $(PKGNAME)
	@rm -r $(PKGNAME)

-include .depend
