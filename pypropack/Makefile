# Makefile for HLSVD
# Created by Brian J. Soher
# Copyright & license info for this Makefile (not HLSVD itself) is in the
# Vespa LICENSE file.

# Find out what platform we are on
ifeq '$(findstring ;,$(PATH))' ';'
    UNAME := Windows
else
    UNAME := $(shell uname 2>/dev/null || echo Unknown)
    UNAME := $(patsubst CYGWIN%,Cygwin,$(UNAME))
    UNAME := $(patsubst MSYS%,MSYS,$(UNAME))
    UNAME := $(patsubst MINGW%,MSYS,$(UNAME))
endif

# If GNU msys is installed, make's default shell will (might?) become sh.exe
# which is Unix-y. That doesn't play well with the DOS-specific commands for
# recreating BIN_DIR.
ifeq ($(UNAME),Windows)
SHELL=cmd
endif

DIST_TARGET_DIR := ../hlsvdpro

PYMAJ := $(shell python -c "from __future__ import print_function;import sys; print(sys.version_info.major)")

all:
ifeq ($(UNAME),Windows)
	@echo bil $(UNAME)
	copy .\complex16\Makefile.win .\complex16\Makefile
	copy .\complex8\Makefile.win .\complex8\Makefile
	copy .\double\Makefile.win .\double\Makefile
	copy .\single\Makefile.win .\single\Makefile
	if EXIST .\complex16\bin rmdir /s /q .\complex16\bin
	if EXIST .\complex8\bin rmdir /s /q .\complex8\bin
	if EXIST .\double\bin rmdir /s /q .\double\bin
	if EXIST .\single\bin rmdir /s /q .\single\bin
	copy .\complex16\f2py_Util\py$(PYMAJ)_f2py\* .\complex16\f2py_Util\.
	copy .\complex8\f2py_Util\py$(PYMAJ)_f2py\* .\complex8\f2py_Util\.
	copy .\double\f2py_Util\py$(PYMAJ)_f2py\* .\double\f2py_Util\.
	copy .\single\f2py_Util\py$(PYMAJ)_f2py\* .\single\f2py_Util\.
	cd .\complex16 && $(MAKE)
	cd .\complex8 && $(MAKE)
	cd .\double && $(MAKE)
	cd .\single && $(MAKE)
endif
ifeq ($(UNAME),Darwin)        # Mac OS X
	@echo dar $(UNAME) 
	cp complex16/Makefile.osx complex16/Makefile
	cp complex8/Makefile.osx complex8/Makefile
	cp double/Makefile.osx double/Makefile
	cp single/Makefile.osx single/Makefile	
	rm -rf complex16/bin
	rm -rf complex8/bin
	rm -rf double/bin
	rm -rf single/bin
	cp ./complex16/f2py_Util/py$(PYMAJ)_f2py/* ./complex16/f2py_Util/.
	cp ./complex8/f2py_Util/py$(PYMAJ)_f2py/* ./complex8/f2py_Util/.
	cp ./double/f2py_Util/py$(PYMAJ)_f2py/* ./double/f2py_Util/.
	cp ./single/f2py_Util/py$(PYMAJ)_f2py/* ./single/f2py_Util/.
	cd ./complex16 && $(MAKE)
	cd ./complex8 && $(MAKE)
	cd ./double && $(MAKE)
	cd ./single && $(MAKE)
endif
ifeq ($(UNAME),Linux)
	@echo lin $(UNAME)
	cp complex16/Makefile.linux complex16/Makefile
	cp complex8/Makefile.linux complex8/Makefile
	cp double/Makefile.linux double/Makefile
	cp single/Makefile.linux single/Makefile	
	rm -rf complex16/bin
	rm -rf complex8/bin
	rm -rf double/bin
	rm -rf single/bin
	cp ./complex16/f2py_Util/py$(PYMAJ)_f2py/* ./complex16/f2py_Util/.
	cp ./complex8/f2py_Util/py$(PYMAJ)_f2py/* ./complex8/f2py_Util/.
	cp ./double/f2py_Util/py$(PYMAJ)_f2py/* ./double/f2py_Util/.
	cp ./single/f2py_Util/py$(PYMAJ)_f2py/* ./single/f2py_Util/.
	cd ./complex16 && $(MAKE)
	cd ./complex8 && $(MAKE)
	cd ./double && $(MAKE)
	cd ./single && $(MAKE)
endif


dist:
ifeq ($(UNAME),Windows)
	@echo bil $(UNAME)
	if EXIST ..\hlsvdpro del /f /q ..\hlsvdpro\*propack.pyd
	copy .\complex16\Makefile.win .\complex16\Makefile
	copy .\complex8\Makefile.win .\complex8\Makefile
	copy .\double\Makefile.win .\double\Makefile
	copy .\single\Makefile.win .\single\Makefile
	cd .\complex16 && $(MAKE) dist
	cd .\complex8 && $(MAKE) dist
	cd .\double && $(MAKE) dist
	cd .\single && $(MAKE) dist
endif
ifeq ($(UNAME),Darwin)        # Mac OS X
	@echo dar $(UNAME) 
	rm -rf $(DIST_TARGET_DIR)/*propack.dylib
	cp complex16/Makefile.osx complex16/Makefile
	cp complex8/Makefile.osx complex8/Makefile
	cp double/Makefile.osx double/Makefile
	cp single/Makefile.osx single/Makefile
	cd ./complex16 && $(MAKE) dist
	cd ./complex8 && $(MAKE) dist
	cd ./double && $(MAKE) dist
	cd ./single && $(MAKE) dist
endif
ifeq ($(UNAME),Linux)
	@echo lin $(UNAME)
	rm -rf $(DIST_TARGET_DIR)/*propack.so
	cp complex16/Makefile.linux complex16/Makefile
	cp complex8/Makefile.linux complex8/Makefile
	cp double/Makefile.linux double/Makefile
	cp single/Makefile.linux single/Makefile
	cd ./complex16 && $(MAKE) dist
	cd ./complex8 && $(MAKE) dist
	cd ./double && $(MAKE) dist
	cd ./single && $(MAKE) dist
endif


clean:
ifeq ($(UNAME),Windows)
	@echo bil $(UNAME)
	copy .\complex16\Makefile.win .\complex16\Makefile
	copy .\complex8\Makefile.win .\complex8\Makefile
	copy .\double\Makefile.win .\double\Makefile
	copy .\single\Makefile.win .\single\Makefile
	cd .\complex16 && $(MAKE) clean
	cd .\complex8 && $(MAKE) clean
	cd .\double && $(MAKE) clean
	cd .\single && $(MAKE) clean
endif
ifeq ($(UNAME),Darwin)        # Mac OS X
	@echo dar $(UNAME) 
	cp complex16/Makefile.osx complex16/Makefile
	cp complex8/Makefile.osx complex8/Makefile
	cp double/Makefile.osx double/Makefile
	cp single/Makefile.osx single/Makefile
	cd ./complex16 && $(MAKE) clean
	cd ./complex8 && $(MAKE) clean
	cd ./double && $(MAKE) clean
	cd ./single && $(MAKE) clean
endif
ifeq ($(UNAME),Linux)
	@echo lin $(UNAME)
	cp complex16/Makefile.linux complex16/Makefile
	cp complex8/Makefile.linux complex8/Makefile
	cp double/Makefile.linux double/Makefile
	cp single/Makefile.linux single/Makefile
	cd ./complex16 && $(MAKE) clean
	cd ./complex8 && $(MAKE) clean
	cd ./double && $(MAKE) clean
	cd ./single && $(MAKE) clean
endif
