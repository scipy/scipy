Name: nifticlib
Summary: Niftilib C library
Version: 0.5
Release: 1
License: Public Domain
Group: Applications/Scientific
Source: %{name}-%{version}.tar.gz
#Patch1: nifticlib.patch
URL: http://niftilib.sourceforge.net/
BuildRoot: %{_tmppath}/%{name}-%{version}-root
Packager: Andy Loening <loening at alum dot mit dot edu>
#Requires: openssl
BuildRequires: cmake


%description
Nifticlib is a set of C i/o libraries for reading and writing files in
the nifti-1 data format. nifti-1 is a binary file format for storing
medical image data, e.g. magnetic resonance image (MRI) and functional
MRI (fMRI) brain images.

%package devel
Summary: static libraries and header files for nifticlib development
Group: Development/Libraries
Requires: nifticlib = %{version}

%description devel

The nifticlib-devel package contains the header files and static libraries
necessary for developing programs that make use of the nifticlib library.

%prep
%setup -q

%build
cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=%{_prefix} -DNIFTI_INSTALL_INCLUDE_DIR=%{_prefix}/include/nifti  .
make

%install
rm -rf $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT

## hack to get this to work for x86_64
%if "%{?_lib}" == "lib64" 
	mv $RPM_BUILD_ROOT/usr/lib $RPM_BUILD_ROOT/%{_libdir}
%endif


%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%doc README
%{_libdir}/*
#%{_datadir}/doc/*

%files devel
%defattr(-,root,root)
%{_bindir}/*
%{_includedir}/*

