#
# Copyright 2018 The Statslabs Authors.
#
cmake_minimum_required(VERSION 3.0)
project(rmath VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)

##################################################
# Create target and set properties
##################################################

add_library(rmath
        src/mlutils.c
        src/d1mach.c src/i1mach.c
        src/fmax2.c src/fmin2.c src/fprec.c src/fround.c src/ftrunc.c
        src/sign.c src/fsign.c
        src/imax2.c src/imin2.c
        src/chebyshev.c src/log1p.c src/expm1.c src/lgammacor.c src/gammalims.c src/stirlerr.c src/bd0.c
        src/gamma.c src/lgamma.c src/gamma_cody.c src/beta.c src/lbeta.c src/polygamma.c src/cospi.c
        src/bessel_i.c src/bessel_j.c src/bessel_k.c src/bessel_y.c
        src/choose.c
        src/snorm.c src/sexp.c
        src/dgamma.c src/pgamma.c src/qgamma.c src/rgamma.c
        src/dbeta.c src/pbeta.c src/qbeta.c src/rbeta.c
        src/dunif.c src/punif.c src/qunif.c src/runif.c
        src/dnorm.c src/pnorm.c src/qnorm.c src/rnorm.c
        src/dlnorm.c src/plnorm.c src/qlnorm.c src/rlnorm.c
        src/df.c src/pf.c src/qf.c src/rf.c src/dnf.c
        src/dt.c src/pt.c src/qt.c src/rt.c src/dnt.c
        src/dchisq.c src/pchisq.c src/qchisq.c src/rchisq.c src/rnchisq.c
        src/dbinom.c src/pbinom.c src/qbinom.c src/rbinom.c
        src/rmultinom.c
        src/dcauchy.c src/pcauchy.c src/qcauchy.c src/rcauchy.c
        src/dexp.c src/pexp.c src/qexp.c src/rexp.c
        src/dgeom.c src/pgeom.c src/qgeom.c src/rgeom.c
        src/dhyper.c src/phyper.c src/qhyper.c src/rhyper.c
        src/dnbinom.c src/pnbinom.c src/qnbinom.c src/rnbinom.c
        src/dpois.c src/ppois.c src/qpois.c src/rpois.c
        src/dweibull.c src/pweibull.c src/qweibull.c src/rweibull.c
        src/dlogis.c src/plogis.c src/qlogis.c src/rlogis.c
        src/dnchisq.c src/pnchisq.c src/qnchisq.c
        src/dnbeta.c src/pnbeta.c src/qnbeta.c
        src/pnf.c src/pnt.c src/qnf.c src/qnt.c
        src/ptukey.c src/qtukey.c
        src/toms708.c src/wilcox.c
        src/signrank.c src/sunif.c)

target_link_libraries(rmath m)

# Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(Rmath::rmath ALIAS rmath)

# Set target properties
target_include_directories(rmath
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

target_compile_definitions(rmath PUBLIC MATHLIB_STANDALONE)

##################################################
# Installation instructions
##################################################

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Rmath)

install(TARGETS rmath
        EXPORT rmath-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

#This is required so that the exported target has the name Rmath and not rmath
set_target_properties(rmath PROPERTIES EXPORT_NAME Rmath)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export the targets to a script
install(EXPORT rmath-targets
        FILE
        RmathTargets.cmake
        NAMESPACE
        Rmath::
        DESTINATION
        ${INSTALL_CONFIGDIR}
        )

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/RmathConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/RmathConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/RmathConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
        )

#Install the config, configversion and custom find modules
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/RmathConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/RmathConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )

export(EXPORT rmath-targets FILE ${CMAKE_CURRENT_BINARY_DIR}/RmathTargets.cmake NAMESPACE Rmath::)

# Register package in user's package registry
export(PACKAGE Rmath)
