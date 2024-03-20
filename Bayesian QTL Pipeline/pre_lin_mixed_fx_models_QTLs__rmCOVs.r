
#### Import data
#qtl_dat <- read.table("/Users/dc/Desktop/rsch/pec2/bayesian_lin_mixed_FX_models/AC107223.1__chr4:143633889:A:G__input_for_mxd_fx_model.dat", header=TRUE)
#qtl_dat <- read.table("/Users/dc/Desktop/rsch/pec2/v5_scqtls/bayesian_lin_mixed_FX_models/OR2T10__chr1:248678345:G:A__input_for_mxd_fx_model.dat", header=TRUE)

args = commandArgs(trailingOnly=TRUE)
input_file_name = args[1]
qtl_dat <- read.table(input_file_name, header=TRUE)


num_covariates = 128
col_index_w_expr_data = num_covariates + 3


ids<-sort(unique(qtl_dat$cellType))
m<-length(ids)
Y<-list()   ;   X<-list()   ;   N<-NULL

COV_1<-list()   ;   COV_2<-list()   ;   COV_3<-list()   ;   COV_4<-list()   ;   COV_5<-list()   ;   COV_6<-list()   ;   COV_7<-list()   ;   COV_8<-list()   ;   COV_9<-list()   ;   COV_10<-list()
COV_11<-list()   ;   COV_12<-list()   ;   COV_13<-list()   ;   COV_14<-list()   ;   COV_15<-list()   ;   COV_16<-list()   ;   COV_17<-list()   ;   COV_18<-list()   ;   COV_19<-list()   ;   COV_20<-list()
COV_21<-list()   ;   COV_22<-list()   ;   COV_23<-list()   ;   COV_24<-list()   ;   COV_25<-list()   ;   COV_26<-list()   ;   COV_27<-list()   ;   COV_28<-list()   ;   COV_29<-list()   ;   COV_30<-list()
COV_31<-list()   ;   COV_32<-list()   ;   COV_33<-list()   ;   COV_34<-list()   ;   COV_35<-list()   ;   COV_36<-list()   ;   COV_37<-list()   ;   COV_38<-list()   ;   COV_39<-list()   ;   COV_40<-list()
COV_41<-list()   ;   COV_42<-list()   ;   COV_43<-list()   ;   COV_44<-list()   ;   COV_45<-list()   ;   COV_46<-list()   ;   COV_47<-list()   ;   COV_48<-list()   ;   COV_49<-list()   ;   COV_50<-list()
COV_51<-list()   ;   COV_52<-list()   ;   COV_53<-list()   ;   COV_54<-list()   ;   COV_55<-list()   ;   COV_56<-list()   ;   COV_57<-list()   ;   COV_58<-list()   ;   COV_59<-list()   ;   COV_60<-list()
COV_61<-list()   ;   COV_62<-list()   ;   COV_63<-list()   ;   COV_64<-list()   ;   COV_65<-list()   ;   COV_66<-list()   ;   COV_67<-list()   ;   COV_68<-list()   ;   COV_69<-list()   ;   COV_70<-list()
COV_71<-list()   ;   COV_72<-list()   ;   COV_73<-list()   ;   COV_74<-list()   ;   COV_75<-list()   ;   COV_76<-list()   ;   COV_77<-list()   ;   COV_78<-list()   ;   COV_79<-list()   ;   COV_80<-list()
COV_81<-list()   ;   COV_82<-list()   ;   COV_83<-list()   ;   COV_84<-list()   ;   COV_85<-list()   ;   COV_86<-list()   ;   COV_87<-list()   ;   COV_88<-list()   ;   COV_89<-list()   ;   COV_90<-list()
COV_91<-list()   ;   COV_92<-list()   ;   COV_93<-list()   ;   COV_94<-list()   ;   COV_95<-list()   ;   COV_96<-list()   ;   COV_97<-list()   ;   COV_98<-list()   ;   COV_99<-list()   ;   COV_100<-list()
COV_101<-list()   ;   COV_102<-list()   ;   COV_103<-list()   ;   COV_104<-list()   ;   COV_105<-list()   ;   COV_106<-list()   ;   COV_107<-list()   ;   COV_108<-list()   ;   COV_109<-list()   ;   COV_110<-list()
COV_111<-list()   ;   COV_112<-list()   ;   COV_113<-list()   ;   COV_114<-list()   ;   COV_115<-list()   ;   COV_116<-list()   ;   COV_117<-list()   ;   COV_118<-list()   ;   COV_119<-list()   ;   COV_120<-list()
COV_121<-list()   ;   COV_122<-list()   ;   COV_123<-list()   ;   COV_124<-list()   ;   COV_125<-list()   ;   COV_126<-list()   ;   COV_127<-list()   ;   COV_128<-list()


for(j in 1:m)
{
	Y[[j]]<-qtl_dat[qtl_dat$cellType==ids[j], col_index_w_expr_data]

	N[j]<- sum(qtl_dat$cellType==ids[j])

	xj<-qtl_dat[qtl_dat$cellType==ids[j], 2]

	X[[j]]<-cbind( rep(1,N[j]), xj )

	cov_1_j<-qtl_dat[qtl_dat$cellType==ids[j], 3]
	COV_1[[j]]<-cbind( cov_1_j )

	cov_2_j<-qtl_dat[qtl_dat$cellType==ids[j], 4]
	COV_2[[j]]<-cbind( cov_2_j )

	cov_3_j<-qtl_dat[qtl_dat$cellType==ids[j], 5]
	COV_3[[j]]<-cbind( cov_3_j )

	cov_4_j<-qtl_dat[qtl_dat$cellType==ids[j], 6]
	COV_4[[j]]<-cbind( cov_4_j )

	cov_5_j<-qtl_dat[qtl_dat$cellType==ids[j], 7]
	COV_5[[j]]<-cbind( cov_5_j )

	cov_6_j<-qtl_dat[qtl_dat$cellType==ids[j], 8]
	COV_6[[j]]<-cbind( cov_6_j )

	cov_7_j<-qtl_dat[qtl_dat$cellType==ids[j], 9]
	COV_7[[j]]<-cbind( cov_7_j )

	cov_8_j<-qtl_dat[qtl_dat$cellType==ids[j], 10]
	COV_8[[j]]<-cbind( cov_8_j )

	cov_9_j<-qtl_dat[qtl_dat$cellType==ids[j], 11]
	COV_9[[j]]<-cbind( cov_9_j )

	cov_10_j<-qtl_dat[qtl_dat$cellType==ids[j], 12]
	COV_10[[j]]<-cbind( cov_10_j )

	cov_11_j<-qtl_dat[qtl_dat$cellType==ids[j], 13]
	COV_11[[j]]<-cbind( cov_11_j )

	cov_12_j<-qtl_dat[qtl_dat$cellType==ids[j], 14]
	COV_12[[j]]<-cbind( cov_12_j )

	cov_13_j<-qtl_dat[qtl_dat$cellType==ids[j], 15]
	COV_13[[j]]<-cbind( cov_13_j )

	cov_14_j<-qtl_dat[qtl_dat$cellType==ids[j], 16]
	COV_14[[j]]<-cbind( cov_14_j )

	cov_15_j<-qtl_dat[qtl_dat$cellType==ids[j], 17]
	COV_15[[j]]<-cbind( cov_15_j )

	cov_16_j<-qtl_dat[qtl_dat$cellType==ids[j], 18]
	COV_16[[j]]<-cbind( cov_16_j )

	cov_17_j<-qtl_dat[qtl_dat$cellType==ids[j], 19]
	COV_17[[j]]<-cbind( cov_17_j )

	cov_18_j<-qtl_dat[qtl_dat$cellType==ids[j], 20]
	COV_18[[j]]<-cbind( cov_18_j )

	cov_19_j<-qtl_dat[qtl_dat$cellType==ids[j], 21]
	COV_19[[j]]<-cbind( cov_19_j )

	cov_20_j<-qtl_dat[qtl_dat$cellType==ids[j], 22]
	COV_20[[j]]<-cbind( cov_20_j )

	cov_21_j<-qtl_dat[qtl_dat$cellType==ids[j], 23]
	COV_21[[j]]<-cbind( cov_21_j )

	cov_22_j<-qtl_dat[qtl_dat$cellType==ids[j], 24]
	COV_22[[j]]<-cbind( cov_22_j )

	cov_23_j<-qtl_dat[qtl_dat$cellType==ids[j], 25]
	COV_23[[j]]<-cbind( cov_23_j )

	cov_24_j<-qtl_dat[qtl_dat$cellType==ids[j], 26]
	COV_24[[j]]<-cbind( cov_24_j )

	cov_25_j<-qtl_dat[qtl_dat$cellType==ids[j], 27]
	COV_25[[j]]<-cbind( cov_25_j )

	cov_26_j<-qtl_dat[qtl_dat$cellType==ids[j], 28]
	COV_26[[j]]<-cbind( cov_26_j )

	cov_27_j<-qtl_dat[qtl_dat$cellType==ids[j], 29]
	COV_27[[j]]<-cbind( cov_27_j )

	cov_28_j<-qtl_dat[qtl_dat$cellType==ids[j], 30]
	COV_28[[j]]<-cbind( cov_28_j )

	cov_29_j<-qtl_dat[qtl_dat$cellType==ids[j], 31]
	COV_29[[j]]<-cbind( cov_29_j )

	cov_30_j<-qtl_dat[qtl_dat$cellType==ids[j], 32]
	COV_30[[j]]<-cbind( cov_30_j )

	cov_31_j<-qtl_dat[qtl_dat$cellType==ids[j], 33]
	COV_31[[j]]<-cbind( cov_31_j )

	cov_32_j<-qtl_dat[qtl_dat$cellType==ids[j], 34]
	COV_32[[j]]<-cbind( cov_32_j )

	cov_33_j<-qtl_dat[qtl_dat$cellType==ids[j], 35]
	COV_33[[j]]<-cbind( cov_33_j )

	cov_34_j<-qtl_dat[qtl_dat$cellType==ids[j], 36]
	COV_34[[j]]<-cbind( cov_34_j )

	cov_35_j<-qtl_dat[qtl_dat$cellType==ids[j], 37]
	COV_35[[j]]<-cbind( cov_35_j )

	cov_36_j<-qtl_dat[qtl_dat$cellType==ids[j], 38]
	COV_36[[j]]<-cbind( cov_36_j )

	cov_37_j<-qtl_dat[qtl_dat$cellType==ids[j], 39]
	COV_37[[j]]<-cbind( cov_37_j )

	cov_38_j<-qtl_dat[qtl_dat$cellType==ids[j], 40]
	COV_38[[j]]<-cbind( cov_38_j )

	cov_39_j<-qtl_dat[qtl_dat$cellType==ids[j], 41]
	COV_39[[j]]<-cbind( cov_39_j )

	cov_40_j<-qtl_dat[qtl_dat$cellType==ids[j], 42]
	COV_40[[j]]<-cbind( cov_40_j )

	cov_41_j<-qtl_dat[qtl_dat$cellType==ids[j], 43]
	COV_41[[j]]<-cbind( cov_41_j )

	cov_42_j<-qtl_dat[qtl_dat$cellType==ids[j], 44]
	COV_42[[j]]<-cbind( cov_42_j )

	cov_43_j<-qtl_dat[qtl_dat$cellType==ids[j], 45]
	COV_43[[j]]<-cbind( cov_43_j )

	cov_44_j<-qtl_dat[qtl_dat$cellType==ids[j], 46]
	COV_44[[j]]<-cbind( cov_44_j )

	cov_45_j<-qtl_dat[qtl_dat$cellType==ids[j], 47]
	COV_45[[j]]<-cbind( cov_45_j )

	cov_46_j<-qtl_dat[qtl_dat$cellType==ids[j], 48]
	COV_46[[j]]<-cbind( cov_46_j )

	cov_47_j<-qtl_dat[qtl_dat$cellType==ids[j], 49]
	COV_47[[j]]<-cbind( cov_47_j )

	cov_48_j<-qtl_dat[qtl_dat$cellType==ids[j], 50]
	COV_48[[j]]<-cbind( cov_48_j )

	cov_49_j<-qtl_dat[qtl_dat$cellType==ids[j], 51]
	COV_49[[j]]<-cbind( cov_49_j )

	cov_50_j<-qtl_dat[qtl_dat$cellType==ids[j], 52]
	COV_50[[j]]<-cbind( cov_50_j )

	cov_51_j<-qtl_dat[qtl_dat$cellType==ids[j], 53]
	COV_51[[j]]<-cbind( cov_51_j )

	cov_52_j<-qtl_dat[qtl_dat$cellType==ids[j], 54]
	COV_52[[j]]<-cbind( cov_52_j )

	cov_53_j<-qtl_dat[qtl_dat$cellType==ids[j], 55]
	COV_53[[j]]<-cbind( cov_53_j )

	cov_54_j<-qtl_dat[qtl_dat$cellType==ids[j], 56]
	COV_54[[j]]<-cbind( cov_54_j )

	cov_55_j<-qtl_dat[qtl_dat$cellType==ids[j], 57]
	COV_55[[j]]<-cbind( cov_55_j )

	cov_56_j<-qtl_dat[qtl_dat$cellType==ids[j], 58]
	COV_56[[j]]<-cbind( cov_56_j )

	cov_57_j<-qtl_dat[qtl_dat$cellType==ids[j], 59]
	COV_57[[j]]<-cbind( cov_57_j )

	cov_58_j<-qtl_dat[qtl_dat$cellType==ids[j], 60]
	COV_58[[j]]<-cbind( cov_58_j )

	cov_59_j<-qtl_dat[qtl_dat$cellType==ids[j], 61]
	COV_59[[j]]<-cbind( cov_59_j )

	cov_60_j<-qtl_dat[qtl_dat$cellType==ids[j], 62]
	COV_60[[j]]<-cbind( cov_60_j )

	cov_61_j<-qtl_dat[qtl_dat$cellType==ids[j], 63]
	COV_61[[j]]<-cbind( cov_61_j )

	cov_62_j<-qtl_dat[qtl_dat$cellType==ids[j], 64]
	COV_62[[j]]<-cbind( cov_62_j )

	cov_63_j<-qtl_dat[qtl_dat$cellType==ids[j], 65]
	COV_63[[j]]<-cbind( cov_63_j )

	cov_64_j<-qtl_dat[qtl_dat$cellType==ids[j], 66]
	COV_64[[j]]<-cbind( cov_64_j )

	cov_65_j<-qtl_dat[qtl_dat$cellType==ids[j], 67]
	COV_65[[j]]<-cbind( cov_65_j )

	cov_66_j<-qtl_dat[qtl_dat$cellType==ids[j], 68]
	COV_66[[j]]<-cbind( cov_66_j )

	cov_67_j<-qtl_dat[qtl_dat$cellType==ids[j], 69]
	COV_67[[j]]<-cbind( cov_67_j )

	cov_68_j<-qtl_dat[qtl_dat$cellType==ids[j], 70]
	COV_68[[j]]<-cbind( cov_68_j )

	cov_69_j<-qtl_dat[qtl_dat$cellType==ids[j], 71]
	COV_69[[j]]<-cbind( cov_69_j )

	cov_70_j<-qtl_dat[qtl_dat$cellType==ids[j], 72]
	COV_70[[j]]<-cbind( cov_70_j )

	cov_71_j<-qtl_dat[qtl_dat$cellType==ids[j], 73]
	COV_71[[j]]<-cbind( cov_71_j )

	cov_72_j<-qtl_dat[qtl_dat$cellType==ids[j], 74]
	COV_72[[j]]<-cbind( cov_72_j )

	cov_73_j<-qtl_dat[qtl_dat$cellType==ids[j], 75]
	COV_73[[j]]<-cbind( cov_73_j )

	cov_74_j<-qtl_dat[qtl_dat$cellType==ids[j], 76]
	COV_74[[j]]<-cbind( cov_74_j )

	cov_75_j<-qtl_dat[qtl_dat$cellType==ids[j], 77]
	COV_75[[j]]<-cbind( cov_75_j )

	cov_76_j<-qtl_dat[qtl_dat$cellType==ids[j], 78]
	COV_76[[j]]<-cbind( cov_76_j )

	cov_77_j<-qtl_dat[qtl_dat$cellType==ids[j], 79]
	COV_77[[j]]<-cbind( cov_77_j )

	cov_78_j<-qtl_dat[qtl_dat$cellType==ids[j], 80]
	COV_78[[j]]<-cbind( cov_78_j )

	cov_79_j<-qtl_dat[qtl_dat$cellType==ids[j], 81]
	COV_79[[j]]<-cbind( cov_79_j )

	cov_80_j<-qtl_dat[qtl_dat$cellType==ids[j], 82]
	COV_80[[j]]<-cbind( cov_80_j )

	cov_81_j<-qtl_dat[qtl_dat$cellType==ids[j], 83]
	COV_81[[j]]<-cbind( cov_81_j )

	cov_82_j<-qtl_dat[qtl_dat$cellType==ids[j], 84]
	COV_82[[j]]<-cbind( cov_82_j )

	cov_83_j<-qtl_dat[qtl_dat$cellType==ids[j], 85]
	COV_83[[j]]<-cbind( cov_83_j )

	cov_84_j<-qtl_dat[qtl_dat$cellType==ids[j], 86]
	COV_84[[j]]<-cbind( cov_84_j )

	cov_85_j<-qtl_dat[qtl_dat$cellType==ids[j], 87]
	COV_85[[j]]<-cbind( cov_85_j )

	cov_86_j<-qtl_dat[qtl_dat$cellType==ids[j], 88]
	COV_86[[j]]<-cbind( cov_86_j )

	cov_87_j<-qtl_dat[qtl_dat$cellType==ids[j], 89]
	COV_87[[j]]<-cbind( cov_87_j )

	cov_88_j<-qtl_dat[qtl_dat$cellType==ids[j], 90]
	COV_88[[j]]<-cbind( cov_88_j )

	cov_89_j<-qtl_dat[qtl_dat$cellType==ids[j], 91]
	COV_89[[j]]<-cbind( cov_89_j )

	cov_90_j<-qtl_dat[qtl_dat$cellType==ids[j], 92]
	COV_90[[j]]<-cbind( cov_90_j )

	cov_91_j<-qtl_dat[qtl_dat$cellType==ids[j], 93]
	COV_91[[j]]<-cbind( cov_91_j )

	cov_92_j<-qtl_dat[qtl_dat$cellType==ids[j], 94]
	COV_92[[j]]<-cbind( cov_92_j )

	cov_93_j<-qtl_dat[qtl_dat$cellType==ids[j], 95]
	COV_93[[j]]<-cbind( cov_93_j )

	cov_94_j<-qtl_dat[qtl_dat$cellType==ids[j], 96]
	COV_94[[j]]<-cbind( cov_94_j )

	cov_95_j<-qtl_dat[qtl_dat$cellType==ids[j], 97]
	COV_95[[j]]<-cbind( cov_95_j )

	cov_96_j<-qtl_dat[qtl_dat$cellType==ids[j], 98]
	COV_96[[j]]<-cbind( cov_96_j )

	cov_97_j<-qtl_dat[qtl_dat$cellType==ids[j], 99]
	COV_97[[j]]<-cbind( cov_97_j )

	cov_98_j<-qtl_dat[qtl_dat$cellType==ids[j], 100]
	COV_98[[j]]<-cbind( cov_98_j )

	cov_99_j<-qtl_dat[qtl_dat$cellType==ids[j], 101]
	COV_99[[j]]<-cbind( cov_99_j )

	cov_100_j<-qtl_dat[qtl_dat$cellType==ids[j], 102]
	COV_100[[j]]<-cbind( cov_100_j )

	cov_101_j<-qtl_dat[qtl_dat$cellType==ids[j], 103]
	COV_101[[j]]<-cbind( cov_101_j )

	cov_102_j<-qtl_dat[qtl_dat$cellType==ids[j], 104]
	COV_102[[j]]<-cbind( cov_102_j )

	cov_103_j<-qtl_dat[qtl_dat$cellType==ids[j], 105]
	COV_103[[j]]<-cbind( cov_103_j )

	cov_104_j<-qtl_dat[qtl_dat$cellType==ids[j], 106]
	COV_104[[j]]<-cbind( cov_104_j )

	cov_105_j<-qtl_dat[qtl_dat$cellType==ids[j], 107]
	COV_105[[j]]<-cbind( cov_105_j )

	cov_106_j<-qtl_dat[qtl_dat$cellType==ids[j], 108]
	COV_106[[j]]<-cbind( cov_106_j )

	cov_107_j<-qtl_dat[qtl_dat$cellType==ids[j], 109]
	COV_107[[j]]<-cbind( cov_107_j )

	cov_108_j<-qtl_dat[qtl_dat$cellType==ids[j], 110]
	COV_108[[j]]<-cbind( cov_108_j )

	cov_109_j<-qtl_dat[qtl_dat$cellType==ids[j], 111]
	COV_109[[j]]<-cbind( cov_109_j )

	cov_110_j<-qtl_dat[qtl_dat$cellType==ids[j], 112]
	COV_110[[j]]<-cbind( cov_110_j )

	cov_111_j<-qtl_dat[qtl_dat$cellType==ids[j], 113]
	COV_111[[j]]<-cbind( cov_111_j )

	cov_112_j<-qtl_dat[qtl_dat$cellType==ids[j], 114]
	COV_112[[j]]<-cbind( cov_112_j )

	cov_113_j<-qtl_dat[qtl_dat$cellType==ids[j], 115]
	COV_113[[j]]<-cbind( cov_113_j )

	cov_114_j<-qtl_dat[qtl_dat$cellType==ids[j], 116]
	COV_114[[j]]<-cbind( cov_114_j )

	cov_115_j<-qtl_dat[qtl_dat$cellType==ids[j], 117]
	COV_115[[j]]<-cbind( cov_115_j )

	cov_116_j<-qtl_dat[qtl_dat$cellType==ids[j], 118]
	COV_116[[j]]<-cbind( cov_116_j )

	cov_117_j<-qtl_dat[qtl_dat$cellType==ids[j], 119]
	COV_117[[j]]<-cbind( cov_117_j )

	cov_118_j<-qtl_dat[qtl_dat$cellType==ids[j], 120]
	COV_118[[j]]<-cbind( cov_118_j )

	cov_119_j<-qtl_dat[qtl_dat$cellType==ids[j], 121]
	COV_119[[j]]<-cbind( cov_119_j )

	cov_120_j<-qtl_dat[qtl_dat$cellType==ids[j], 122]
	COV_120[[j]]<-cbind( cov_120_j )

	cov_121_j<-qtl_dat[qtl_dat$cellType==ids[j], 123]
	COV_121[[j]]<-cbind( cov_121_j )

	cov_122_j<-qtl_dat[qtl_dat$cellType==ids[j], 124]
	COV_122[[j]]<-cbind( cov_122_j )

	cov_123_j<-qtl_dat[qtl_dat$cellType==ids[j], 125]
	COV_123[[j]]<-cbind( cov_123_j )

	cov_124_j<-qtl_dat[qtl_dat$cellType==ids[j], 126]
	COV_124[[j]]<-cbind( cov_124_j )

	cov_125_j<-qtl_dat[qtl_dat$cellType==ids[j], 127]
	COV_125[[j]]<-cbind( cov_125_j )

	cov_126_j<-qtl_dat[qtl_dat$cellType==ids[j], 128]
	COV_126[[j]]<-cbind( cov_126_j )

	cov_127_j<-qtl_dat[qtl_dat$cellType==ids[j], 129]
	COV_127[[j]]<-cbind( cov_127_j )

	cov_128_j<-qtl_dat[qtl_dat$cellType==ids[j], 130]
	COV_128[[j]]<-cbind( cov_128_j )
}



list_of_COV_indeces_to_exclude <- list()
for(j in 1:m)
{
	fit<-lm(Y[[j]] ~ -1 + X[[j]] + COV_1[[j]] + COV_2[[j]] + COV_3[[j]] + COV_4[[j]] + COV_5[[j]] + COV_6[[j]] + COV_7[[j]] + COV_8[[j]] + COV_9[[j]] + COV_10[[j]] + COV_11[[j]] + COV_12[[j]] + COV_13[[j]] + COV_14[[j]] + COV_15[[j]] + COV_16[[j]] + COV_17[[j]] + COV_18[[j]] + COV_19[[j]] + COV_20[[j]] + COV_21[[j]] + COV_22[[j]] + COV_23[[j]] + COV_24[[j]] + COV_25[[j]] + COV_26[[j]] + COV_27[[j]] + COV_28[[j]] + COV_29[[j]] + COV_30[[j]] + COV_31[[j]] + COV_32[[j]] + COV_33[[j]] + COV_34[[j]] + COV_35[[j]] + COV_36[[j]] + COV_37[[j]] + COV_38[[j]] + COV_39[[j]] + COV_40[[j]] + COV_41[[j]] + COV_42[[j]] + COV_43[[j]] + COV_44[[j]] + COV_45[[j]] + COV_46[[j]] + COV_47[[j]] + COV_48[[j]] + COV_49[[j]] + COV_50[[j]] + COV_51[[j]] + COV_52[[j]] + COV_53[[j]] + COV_54[[j]] + COV_55[[j]] + COV_56[[j]] + COV_57[[j]] + COV_58[[j]] + COV_59[[j]] + COV_60[[j]] + COV_61[[j]] + COV_62[[j]] + COV_63[[j]] + COV_64[[j]] + COV_65[[j]] + COV_66[[j]] + COV_67[[j]] + COV_68[[j]] + COV_69[[j]] + COV_70[[j]] + COV_71[[j]] + COV_72[[j]] + COV_73[[j]] + COV_74[[j]] + COV_75[[j]] + COV_76[[j]] + COV_77[[j]] + COV_78[[j]] + COV_79[[j]] + COV_80[[j]] + COV_81[[j]] + COV_82[[j]] + COV_83[[j]] + COV_84[[j]] + COV_85[[j]] + COV_86[[j]] + COV_87[[j]] + COV_88[[j]] + COV_89[[j]] + COV_90[[j]] + COV_91[[j]] + COV_92[[j]] + COV_93[[j]] + COV_94[[j]] + COV_95[[j]] + COV_96[[j]] + COV_97[[j]] + COV_98[[j]] + COV_99[[j]] + COV_100[[j]] + COV_101[[j]] + COV_102[[j]] + COV_103[[j]] + COV_104[[j]] + COV_105[[j]] + COV_106[[j]] + COV_107[[j]] + COV_108[[j]] + COV_109[[j]] + COV_110[[j]] + COV_111[[j]] + COV_112[[j]] + COV_113[[j]] + COV_114[[j]] + COV_115[[j]] + COV_116[[j]] + COV_117[[j]] + COV_118[[j]] + COV_119[[j]] + COV_120[[j]] + COV_121[[j]] + COV_122[[j]] + COV_123[[j]] + COV_124[[j]] + COV_125[[j]] + COV_126[[j]] + COV_127[[j]] + COV_128[[j]])
	fit$coefficients
	for(c in 3:(num_covariates+2)) {
		if (is.na(fit$coefficients[c])) {
			list_of_COV_indeces_to_exclude <- append(list_of_COV_indeces_to_exclude, c)
		}
	}
}


#list_of_COV_indeces_to_retain <- list()
for ( c in 3:(num_covariates+2) ) {
	if ( !(c %in% list_of_COV_indeces_to_exclude) ){
		#list_of_COV_indeces_to_retain <- append(list_of_COV_indeces_to_retain, c)
		cov_index = c - 2
		cat(  paste0(cov_index, "\n")  )
	}
}

