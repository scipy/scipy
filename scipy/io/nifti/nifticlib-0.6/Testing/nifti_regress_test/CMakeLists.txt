# Do all regression tests
PROJECT(REGRESSION_TEST)

SET(NT ${NIFTI_BINARY_DIR}/utils/nifti_tool )
SET(DATA ${REGRESSION_TEST_BINARY_DIR}/nifti_regress_data)
SET(TESTS ${REGRESSION_TEST_SOURCE_DIR}/cmake_testscripts)

ADD_TEST(fetch_data
sh ${TESTS}/fetch_data_test.sh ${REGRESSION_TEST_BINARY_DIR} 
)

ADD_TEST(nifti_ver ${NT} -nifti_ver)
ADD_TEST(nifti_tool_ver ${NT} -ver)
ADD_TEST(nifti_help ${NT} -help)
ADD_TEST(nifti_nifti_hist ${NT} -nifti_hist)
ADD_TEST(nifti_hist ${NT} -hist)

ADD_TEST(nifti_disp_hdr ${NT} -disp_hdr -infiles ${DATA}/anat0.nii )

ADD_TEST(nifti_disp_nim ${NT} -disp_nim -infiles ${DATA}/anat0.nii )

ADD_TEST(nifti_disp_ext ${NT} -disp_ext -infiles ${DATA}/anat0.nii )

ADD_TEST(nifti_modhdr_exts
sh ${TESTS}/mod_header_test.sh ${NT} ${DATA} )

ADD_TEST(nifti_bricks_test
sh ${TESTS}/bricks_test.sh ${NT} ${DATA} )

ADD_TEST(nifti_dts_test
sh ${TESTS}/dts_test.sh ${NT} ${DATA} )

ADD_TEST(nifti_dci_test
sh ${TESTS}/dci_test.sh ${NT} ${DATA} )

ADD_TEST(nifti_comment_test
sh ${TESTS}/comment_test.sh ${NT} ${DATA} )
  
ADD_TEST(nifti_header_check ${NT} -check_hdr -infiles ${DATA}/anat0.nii )

ADD_TEST(nifti_nim_check ${NT} -check_nim -infiles ${DATA}/anat0.nii )

ADD_TEST(nifti_dsets_test
sh ${TESTS}/dsets_test.sh ${NT} ${DATA} )

ADD_TEST(nifti_newfiles_test
sh ${TESTS}/newfiles_test.sh ${NT} ${DATA} )
