
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
178
179
180
181
182
183
184
185
186
187
188
189
190
191
192
193
194
195
196
197
198
199
200
201
202
203
204
205
206
207
208
209
210
211
212
213
214
215
216
217
218
219
220
221
222
223
224
225
226
227
228
229
230
231
232
233
234
235
236
237
238
239
240
241
242
243
244
245
246
247
248
249
250
251
252
253
254
255
256
257
258
259
260
261
262
263
264
265
266
267
268
269
270
271
272
273
274
275
276
277
278
279
280
281
282
283
284
285
286
287
288
289
290
291
292
293
294
295
296
297
298
299
300
301
302
303
304
305
306
307
308
309
310
311
312
313
314
315
316
317
318
319
320
321
322
323
324
325
326
327
328
329
330
331
332
333
334
335
336
337
338
339
340
341
342
343
344
345
346
347
348
349
350
351
352
353
354
355
356
357
358
359
360
361
362
363
364
365
366
367
368
369
370
371
372
373
374
375
376
377
378
379
380
381
382
383
384
385
386
387
388
389
390
391
392
393
394
395
396
397
398
399
400
401
402
403
404
405
406
407
408
409
410
411
412
413
414
415
416
417
418
419
420
421
422
423
424
425
426
427
428
429
430
431
432
433
434
435
436
437
438
439
440
441
442
443
444
445
446
447
448
449
450
451
452
453
454
455
456
457
458
459
460
461
462
463
464
465
466
467
468
469
470
471
472
473
474
475
476
477
478
479
480
481
482
483
484
485
486
487
488
489
490
491
492
493
494
495
496
497
498
499
500
501
502
503
504
505
506
507
508
509
510
511
512
513
514
515
516
517
518
519
520
521
522
523
524
525
526
527
528
529
530
531
532
533
534
535
536
537
538
539
540
541
542
543
544
545
546
547
548
549
550
551
552
553
554
555
556
557
558
559
560
561
562
563
564
565
566
567
568
569
570
571
572
573
574
575
576
577
578
579
580
581
582
583
584
585
586
587
588
589
590
591
592
593
594
595
596
597
598
599
600
601
602
603
604
605
606
607
608
609
610
611
612
613
614
615
616
617
618
619
620
621
622
623
624
625
626
627
628
629
630
631
632
633
634
635
636
637
638
639
640
641
642
643
644
645
646
647
648
649
650
651
652
653
654
655
656
657
658
659
660
661
662
663
664
665
666
667
668
669
670
671
672
673
674
675
676
677
678
679
680
681
682
683
684
685
686
687
688
689
690
691
692
693
694
695
696
697
698
699
700
701
702
703
704
705
706
707
708
709
710
711
712
713
714
715
716
717
718
719
720
721
722
723
724
725
726
727
728
729
730
731
732
733
734
735
736
737
738
739
740
741
742
743
744
745
746
747
748
749
750
751
752
753
754
755
756
757
758
759
760
761
762
763
764
765
766
767
768
769
770
771
772
773
774
775
776
777
778
779
780
781
782
783
784
785
786
787
788
789
790
791
792
793
794
795
796
797
798
799
800
801
802
803
804
805
806
807
808
809
810
811
812
813
814
815
816
817
818
819
820
821
822
823
824
825
826
827
828
829
830
831
832
833
834
835
836
837
838
839
840
841
842
843
844
845
846
847
848
849
850
851
852
853
854
855
856
857
858
859
860
861
862
863
864
865
866
867
868
869
870
871
872
873
874
875
876
877
878
879
880
881
882
883
884
885
886
887
888
889
890
891
892
893
894
895
896
897
898
899
900
901
902
903
904
905
906
907
908
909
910
911
912
913
914
915
916
917
918
919
920
921
922
923
924
925
926
927
928
929
930
931
932
933
934
935
936
937
938
939
940
941
942
943
944
945
946
947
948
949
950
951
952
953
954
955
956
957
958
959
960
961
962
963
964
965
966
967
968
969
970
971
972
973
974
975
976
977
978
979
980
981
982
983
984
985
986
987
988
989
990
991
992
993
994
995
996
997
998
999
1000
1001
1002
1003
1004
1005
1006
1007
1008
1009
1010
1011
1012
1013
1014
1015
1016
1017
1018
1019
1020
1021
1022
1023
1024
1025
1026
1027
1028
1029
1030
1031
1032
1033
1034
1035
1036
1037
1038
1039
1040
1041
1042
1043
1044
1045
1046
1047
1048
1049
1050
1051
1052
1053
1054
1055
1056
1057
1058
1059
1060
1061
1062
1063
1064
1065
1066
1067
1068
1069
1070
1071
1072
1073
1074
1075
1076
1077
1078
1079
1080
1081
1082
1083
1084
1085
1086
1087
1088
1089
1090
1091
1092
1093
1094
1095
1096
1097
1098
1099
1100
1101
1102
1103
1104
1105
1106
1107
1108
1109
1110
1111
1112
1113
1114
1115
1116
1117
1118
1119
1120
1121
1122
1123
1124
1125
1126
1127
1128
1129
1130
1131
1132
1133
1134
1135
1136
1137
1138
1139
1140
1141
1142
1143
1144
1145
1146
1147
1148
1149
1150
1151
1152
1153
1154
1155
1156
1157
1158
1159
1160
1161
1162
1163
1164
1165
1166
1167
1168
1169
1170
1171
1172
1173
1174
1175
1176
1177
1178
1179
1180
1181
1182
1183
1184
1185
1186
1187
1188
1189
1190
1191
1192
1193
1194
1195
1196
1197
1198
1199
1200
1201
1202
1203
1204
1205
1206
1207
1208
1209
1210
1211
1212
1213
1214
1215
1216
1217
1218
1219
1220
1221
1222
1223
1224
1225
1226
1227
1228
1229
1230
1231
1232
1233
1234
1235
1236
1237
1238
1239
1240
# 0.0.2 Aug 2019 make compatible with variable taxonomy and mapper names
rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, t, column_name, extra_cols, extra_cols2, i, row_list)

make_ancom_plots_0.0.2 <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level="X.OTU.ID", name_levels=taxonomic_level,var1,var2, correction_level=3,newcol="newcol", the_id = "X.SampleID") {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T) %>% dplyr::select_(.dots = list(the_id,var1))
  
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  
  ## OTU - by inocula
  col1drop <- paste0("-",the_id,".x"); col1drop
  col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("sample2")) %>% mutate(Group=as.factor(as.character(get(var1)))) %>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  
  ancom.OTU <- ANCOM( phyl4, sig = 0.05,multcorr = correction_level, repeated=FALSE )
  
  detected_taxa <- data.frame(OTU=unlist(sapply(X = ancom.OTU$detected, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),x))))%>% filter(OTU!="V1")
  detected_taxa$OTU <- as.character(detected_taxa$OTU)
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa$OTU) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  otu5
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu3[,paste0(taxonomic_level)]
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% inner_join(detected_taxa, by = structure(names=taxonomic_level, "OTU")) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  otu3 <- otu4
  
  rm(i)
  ## run the filtered taxa to make plots
  for (i in 1:dim(otu3)[1]) {
    
    rm(phyl5, column_name)
    column_name <- otu3[,taxonomic_level][i]
    try(phyl5 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl5")==F) {
      try(phyl5 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl5 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl5 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl5, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl5$mean-phyl5$sem)*.9),max(phyl5$mean+phyl5$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=otu3$newcol[i])) 
  }
  
  ## make a list of the plots
  if(dim(detected_taxa)[1]>1) {
    plist <- list(p_1)
    for(q in 2:dim(detected_taxa)[1]) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (dim(detected_taxa)[1]>0) {
    plist <- list(p_1)
  }
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}

make_ancom_covar_plots_0.0.2 <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3 ,correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  for (i in 1:length(rownames(phyl2))) {
    if (rownames(phyl2)[i] == "" ) {
      print(i)
      rownames(phyl2)[i] <- paste0("o__unassigned",i)
    }
  }
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  
  ## OTU - by inocula
  # col1drop <- paste0("-",the_id,".x"); col1drop
  # col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  phyl5[1,]
  map3[1,]
  # comparison_test=ANCOM.main(OTUdat=phyl5,
  #                            Vardat=map3,
  #                            adjusted=T,
  #                            repeated=F,
  #                            main.var="genotype2",
  #                            adj.formula="time",
  #                            repeat.var=NULL,
  #                            longitudinal = FALSE,
  #                            multcorr=3,
  #                            random.formula="~1|genotype2/cage",
  #                            sig=0.05,
  #                            prev.cut=.5)
  # 
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  ## run the filtered taxa to make plots
  # i=1
  # for (i in 1:dim(otu3)[1]) {
  # 	
  # 	rm(phyl6, column_name)
  # 	column_name <- otu3[,taxonomic_level][i]
  # 	try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
  # 	if(exists("phyl6")==F) {
  # 		try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
  # 	}
  # 	
  # 	phyl6
  # 	rm(t)
  # 	for (t in 1:length(extra_cols2)) {
  # 		otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
  # 		otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
  # 	}
  # 	
  # 	assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
  # 				 	geom_bar(position=position_dodge(), stat="identity") +
  # 				 	geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
  # 				 								width=.2,                    # Width of the error bars
  # 				 								position=position_dodge(.9)) +
  # 				 	coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
  # 				 	theme(axis.text=element_text(size=14), 
  # 				 				panel.background = element_blank(),
  # 				 				axis.line = element_line(), 
  # 				 				axis.ticks=element_line(), 
  # 				 				axis.title=element_text(size=16),
  # 				 				title=element_text(size=13)) +
  # 				 	labs(y="relative abundance",x=var1,title=otu3$newcol[i])) 
  # }
  # 
  # 
  
  
  otu3[1,]
  
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}

make_ancom_covar_plots_0.0.2_date <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3 ,correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  for (i in 1:length(rownames(phyl2))) {
    if (rownames(phyl2)[i] == "" ) {
      print(i)
      rownames(phyl2)[i] <- paste0("o__unassigned",i)
    }
  }
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  
  ## OTU - by inocula
  # col1drop <- paste0("-",the_id,".x"); col1drop
  # col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"
  
  map3[1,]
  map3 <- map3 %>% mutate(time = ifelse(time=="base",0,ifelse(time=="T1",1,ifelse(time=="T2",2,ifelse(time=="T3",3,4))))) 
  try(map3 <- map3 %>% mutate(year2 = gsub(pattern = "y", replacement = "", x = year)),T)
  try(map3$year <- as.Date(as.character(map3$year2), "%Y"),T)
  map3[1,]
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  phyl5[1,]
  map3[1,]
  # comparison_test=ANCOM.main(OTUdat=phyl5,
  #                            Vardat=map3,
  #                            adjusted=T,
  #                            repeated=F,
  #                            main.var="genotype2",
  #                            adj.formula="time",
  #                            repeat.var=NULL,
  #                            longitudinal = FALSE,
  #                            multcorr=3,
  #                            random.formula="~1|genotype2/cage",
  #                            sig=0.05,
  #                            prev.cut=.5)
  # 
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  ## run the filtered taxa to make plots
  # i=1
  # for (i in 1:dim(otu3)[1]) {
  # 	
  # 	rm(phyl6, column_name)
  # 	column_name <- otu3[,taxonomic_level][i]
  # 	try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
  # 	if(exists("phyl6")==F) {
  # 		try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
  # 	}
  # 	
  # 	phyl6
  # 	rm(t)
  # 	for (t in 1:length(extra_cols2)) {
  # 		otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
  # 		otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
  # 	}
  # 	
  # 	assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
  # 				 	geom_bar(position=position_dodge(), stat="identity") +
  # 				 	geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
  # 				 								width=.2,                    # Width of the error bars
  # 				 								position=position_dodge(.9)) +
  # 				 	coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
  # 				 	theme(axis.text=element_text(size=14), 
  # 				 				panel.background = element_blank(),
  # 				 				axis.line = element_line(), 
  # 				 				axis.ticks=element_line(), 
  # 				 				axis.title=element_text(size=16),
  # 				 				title=element_text(size=13)) +
  # 				 	labs(y="relative abundance",x=var1,title=otu3$newcol[i])) 
  # }
  # 
  # 
  
  
  otu3[1,]
  
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}

make_ancom_covar_plots_0.0.2_dateCOVAR <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3 ,correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  for (i in 1:length(rownames(phyl2))) {
    if (rownames(phyl2)[i] == "" ) {
      print(i)
      rownames(phyl2)[i] <- paste0("o__unassigned",i)
    }
  }
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$Sample.ID <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$Sample.ID <- gsub("_","", phyl3$Sample.ID)
  map2$Sample.ID <- gsub("_","",map2$X.SampleID)
  map2$Sample.ID <- gsub("-","",map2$Sample.ID)
  
  ## OTU - by inocula
  # col1drop <- paste0("-",the_id,".x"); col1drop
  # col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl5 <- phyl3 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(Sample.ID, everything()) %>% dplyr::select(-X.SampleID) %>% filter(Sample.ID %in% phyl5$Sample.ID)
  colnames(map3)[1] <- "Sample.ID"
  
  map3[1,]
  map3 <- map3 %>% mutate(time = ifelse(time=="base",0,ifelse(time=="t1",1,ifelse(time=="t2",2,ifelse(time=="t3",3,4))))) 
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  
  if (taxonomic_level == "X.OTU.ID") {
    colnames(phyl5) <- paste0("X",colnames(phyl5))
    colnames(phyl5)[1] <- "Sample.ID"
  }
  
  phyl5[1,]
  map3[1,]
  
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  
  
  
  otu3[1,]
  
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}

make_ancom_covar_plots_4var <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3, var4=var4, correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,var4, Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  
  ## OTU - by inocula
  # col1drop <- paste0("-",the_id,".x"); col1drop
  # col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl5 <- phyl3 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("sample2")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}

make_ancom_covar_plots_4var_periodname <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3, var4=var4, correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,var4, Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  rownames(phyl2)[1] <- "o__other"
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  colnames(phyl3) <- rownames(phyl2)
  phyl3 <- phyl3 %>% mutate(X.SampleID=rownames(phyl3))
  phyl3
  phyl3$sample2 <- gsub("X","", phyl3$X.SampleID)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  
  phyl5 <- phyl3 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(map3)[1] <- "Sample.ID"
  phyl5[1]
  map3[1,]
  map3$date2 <- as.Date(as.character(map3$date2), "%m/%d/%y")
  
  phyl4 <- phyl5 %>% inner_join(map2, by=c("Sample.ID")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}


make_ancom_covar_plots_periodname_0.0.2 <- function (file_path, map_path="", mapper_file=mapper_file, taxonomic_level=taxonomic_level,var1=var1,var2=var2,var3=var3 ,correction_level=3,newcol="newcol", the_id = "X.SampleID", adjusted = adjusted, repeated = repeated, main.var=main.var, adj.formula=adj.formula, repeat.var = repeat.var, longitudinal = latitudinal, multcorr=multcorr, sig=sig, prev.cut=prev.cut, Group = Group, name_levels = name_levels, random.formula=NULL) {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv(paste0('taxonomy',map_path,'/taxonomy_forR.csv')), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table(mapper_file,comment.char = "", header=T, fill=T,sep="\t") %>% dplyr::select_(.dots = list(the_id,var1,var2,var3,Group))
  map2[1,]
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  phyl3$sample2 <- gsub("X","", phyl3$X.SampleID)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  
  phyl5 <- phyl3 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(phyl5)[1] <- "Sample.ID"
  
  map3 <- map2 %>% dplyr::select(sample2, everything()) %>% dplyr::select(-X.SampleID)
  colnames(map3)[1] <- "Sample.ID"
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("sample2")) %>% mutate(Group = get(Group)) %>% droplevels()#%>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  rm(comparison_test)
  comparison_test=ANCOM.main(OTUdat=phyl5,
                             Vardat=map3,
                             adjusted=adjusted,
                             repeated=repeated,
                             main.var=main.var,
                             adj.formula=adj.formula,
                             repeat.var=repeat.var,
                             longitudinal = longitudinal,
                             multcorr=multcorr,
                             random.formula=random.formula,
                             sig=sig,
                             prev.cut=prev.cut)
  
  write.csv(comparison_test$W.taxa, paste("ancomexcel_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  dt1 <- data.frame(comparison_test$W.taxa) %>% filter(detected_0.9==T) %>% dplyr::select(otu.names)
  detected_taxa <- data.frame(OTU=unlist(sapply(X = dt1$otu.names, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),substring(x,1)))))	%>% unlist() %>% unname() %>% as.character()
  
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% filter(get(taxonomic_level) %in% detected_taxa) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  
  ## run the filtered taxa to make plots
  # i=1
  # for (i in 1:dim(otu3)[1]) {
  # 	
  # 	rm(phyl6, column_name)
  # 	column_name <- otu3[,taxonomic_level][i]
  # 	try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
  # 	if(exists("phyl6")==F) {
  # 		try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
  # 	}
  # 	
  # 	phyl6
  # 	rm(t)
  # 	for (t in 1:length(extra_cols2)) {
  # 		otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
  # 		otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
  # 	}
  # 	
  # 	assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
  # 				 	geom_bar(position=position_dodge(), stat="identity") +
  # 				 	geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
  # 				 								width=.2,                    # Width of the error bars
  # 				 								position=position_dodge(.9)) +
  # 				 	coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
  # 				 	theme(axis.text=element_text(size=14), 
  # 				 				panel.background = element_blank(),
  # 				 				axis.line = element_line(), 
  # 				 				axis.ticks=element_line(), 
  # 				 				axis.title=element_text(size=16),
  # 				 				title=element_text(size=13)) +
  # 				 	labs(y="relative abundance",x=var1,title=otu3$newcol[i])) 
  # }
  # 
  # 
  
  
  otu3[1,]
  
  
  for (i in 1:length(detected_taxa)) {
    
    rm(phyl6, column_name)
    column_name <- detected_taxa[i]
    try(phyl6 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl6")==F) {
      try(phyl6 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    phyl6
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl6 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl6, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl6$mean-phyl6$sem)*.9),max(phyl6$mean+phyl6$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=	otu3 %>% filter(get(taxonomic_level)== paste0(column_name)) %>% dplyr::select(newcol) %>% unlist() %>% unname()))
  }
  
  p_1
  ## make a list of the plots
  if(length(detected_taxa)>1) {
    plist <- list(p_1)
    for(q in 2:length(detected_taxa)) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (length(detected_taxa)>0) {
    plist <- list(p_1)
  }
  
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}



make_ancom_plots_0.0.1_period_name <- function (file_path, taxonomic_level="X.OTU.ID", name_levels=taxonomic_level,var1, correction_level=3,newcol="newcol", the_id = "X.SampleID") {
  
  rm(otu_table, ref_names, otu2, map2, phyl2, phyl3, col1drop, col2drop, phyl4, ancom.OTU, detected_taxa, otu3, otu4, otu5, phyl5, plist, row_list)
  otu_table <- read.table(paste('core-metrics-results-',file_path,'/rarefied_table.txt',sep=""), comment.char="", header=T, sep="\t", fill=T, skip=1) %>% left_join(read.csv('taxonomy/taxonomy_forR.csv'), by=c("X.OTU.ID"="Feature.ID"))
  #otu_table$OTU <- as.character(otu_table$X.OTU.ID)
  ref_names <- subset(colnames(otu_table), (colnames(otu_table)%in%c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")==F))
  otu2 <- otu_table[,c("X.OTU.ID","kingdom","phylum","class","order","family","genus","species")]
  otu2$OTU <- as.character(otu2$X.OTU.ID)
  
  map2 <- read.table("metadata.tsv",comment.char = "", header=T) %>% dplyr::select_(.dots = list(the_id,var1))
  
  ## hardcoded
  phyl2 <- otu_table %>% group_by_(.dots=taxonomic_level) %>% summarize_at(ref_names,sum, na.rm=T)
  row_list <- as.character(unlist(phyl2[,paste(taxonomic_level)])) 
  rownames(phyl2) <- gsub(" ","", row_list)
  phyl3 <- phyl2 %>% dplyr::select_(.dots=paste0("-",taxonomic_level)) %>% t() %>% data.frame() 
  phyl3 <- mutate(phyl3, X.SampleID=rownames(phyl3))
  #	phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("X","", phyl3$X.SampleID)
  map2$sample2 <- gsub("_","",map2$X.SampleID)
  map2$sample2 <- gsub("-","",map2$sample2)
  
  ## OTU - by inocula
  col1drop <- paste0("-",the_id,".x"); col1drop
  col2drop <- paste0("-",the_id,".y"); col2drop
  
  phyl4 <- phyl3 %>% inner_join(map2, by=c("sample2")) %>% mutate(Group=as.factor(as.character(get(var1)))) %>% dplyr::select_(.dots = list(col1drop, col2drop, paste0("-",var1))) %>% dplyr::select(-sample2)
  
  ancom.OTU <- ANCOM( phyl4, sig = 0.05,multcorr = correction_level, repeated=FALSE )
  
  detected_taxa <- data.frame(OTU=unlist(sapply(X = ancom.OTU$detected, function(x) ifelse(substring(x,1,1)=="X",substring(x,2),x))))
  detected_taxa$OTU <- as.character(detected_taxa$OTU)
  print(detected_taxa$OTU)
  if(taxonomic_level!="OTU") {
    otu2[,paste(taxonomic_level)] = gsub(" ","",otu2[,paste(taxonomic_level)])
  }
  
  ## filter to OTU table by the data frame
  otu3 <- otu2 %>% mutate(newcol = get(name_levels[1]))#%>% dplyr::select_(.dots=c(taxonomic_level)) %>% inner_join(detected_taxa)
  
  ## melt OTU table for non-OTU stuff
  # melted_table <- data.frame(OTU_ID=character(),Count=integer(), Sample=character(), stringsAsFactors = F)
  # for (i in 2:(dim(otu_table)[2]-7-1)) {
  # 	rm(new_table)
  # 	new_table <- data.frame(OTU=as.character(otu_table$OTU),Count=as.integer(otu_table[,i]),Sample=as.character(colnames(otu_table)[i]), stringsAsFactors = F)
  # 	melted_table <- rbind(melted_table, new_table)
  # }
  # 
  # ## start merging
  # melted_table_no0 <- melted_table %>% filter(Count>-1)
  # mt0 <- melted_table_no0 %>% inner_join(otu2, by=c("OTU"="X.OTU.ID"))
  # mt0[1,]
  # mt0$sample2 <- gsub("\\.","", mt0$Sample)
  # mt0$sample2 <- gsub("_","", mt0$sample2)
  # map2$sample2 <- gsub("_","",map2$X.SampleID)
  # map2$sample2 <- gsub("-","",map2$sample2)
  # mt0[1,]
  # mt1 <- mt0 %>% inner_join(map2)
  # mt2 <- mt1 %>% arrange(phylum, class,order,family,genus,species,OTU)
  # mt2$sort_order <- c(1:dim(mt2)[1])
  # mt2$class <- paste(mt2$phylum, mt2$class,sep="_")
  # mt2$order <- paste(mt2$class, mt2$order,sep="_")
  # mt2$family <- paste(mt2$order, mt2$family,sep="_")
  # mt2$genus <- paste(mt2$family, mt2$genus,sep="_")
  # mt2$species <- paste(mt2$genus, mt2$species,sep="_")
  # mt2$OTU <- paste(mt2$species, mt2$OTU,sep="_")
  # 
  otu5 <- otu2 %>% filter(get(taxonomic_level)%in%detected_taxa$OTU) %>% distinct(get(taxonomic_level), .keep_all = T)
  
  extra_cols2 <- names(table(list(phyl4$Group)))
  extra_cols <- unname(c(extra_cols2,sapply(extra_cols2, function(x) paste0(x,"_sem"))))
  otu5[extra_cols] <- "NA"
  
  otu5
  if(length(name_levels)>1) {
    for(i in 2:length(name_levels)) {
      otu3 <- otu3 %>% mutate(newcol = paste(newcol, get(name_levels[i]), sep="_"))
    }
  }
  
  phyl3$sample2 <- gsub("\\.","", phyl3$X.SampleID)
  phyl3$sample2 <- gsub("_","", phyl3$sample2)
  
  otu3[,paste0(taxonomic_level)]
  otu4 <- otu3 %>% mutate(!!taxonomic_level:=gsub("\\[",".", get(taxonomic_level))) %>% mutate(!!taxonomic_level:=gsub("\\]",".", get(taxonomic_level))) %>% dplyr::select_(.dots=c(taxonomic_level, newcol)) %>% inner_join(detected_taxa, by = structure(names=taxonomic_level, "OTU")) %>% distinct_(.dots=taxonomic_level, .keep_all=T)
  otu3 <- otu4
  
  rm(i)
  ## run the filtered taxa to make plots
  for (i in 1:dim(otu3)[1]) {
    
    rm(phyl5, column_name)
    column_name <- otu3[,taxonomic_level][i]
    try(phyl5 <-phyl4 %>% mutate(rabun = get(column_name)/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))),T)
    if(exists("phyl5")==F) {
      try(phyl5 <- phyl4 %>% mutate(rabun = get(paste("X",column_name, sep=""))/sum(otu_table[,2])) %>% dplyr::select(rabun, Group) %>% group_by(Group) %>% summarize(mean=mean(rabun), sem=sd(rabun)/sqrt(length(rabun))))
    }
    
    rm(t)
    for (t in 1:length(extra_cols2)) {
      otu5[i,paste0(extra_cols2[t])] <- phyl5 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(mean) %>% unlist()
      otu5[i,paste0(extra_cols2[t],"_sem")] <- phyl5 %>% filter(Group==extra_cols2[t]) %>% dplyr::select(sem) %>% unlist()
    }
    
    assign(paste("p",i,sep="_"),value = ggplot(phyl5, aes(x=Group, y=mean)) + 
             geom_bar(position=position_dodge(), stat="identity") +
             geom_errorbar(aes(ymin=mean+sem, ymax=mean-sem),
                           width=.2,                    # Width of the error bars
                           position=position_dodge(.9)) +
             coord_cartesian(ylim=c((min(phyl5$mean-phyl5$sem)*.9),max(phyl5$mean+phyl5$sem)*1.3)) + #scale_y_continuous(limits=c(.25,.4)) + 
             theme(axis.text=element_text(size=14), 
                   panel.background = element_blank(),
                   axis.line = element_line(), 
                   axis.ticks=element_line(), 
                   axis.title=element_text(size=16),
                   title=element_text(size=13)) +
             labs(y="relative abundance",x=var1,title=otu3$newcol[i])) 
  }
  
  ## make a list of the plots
  if(dim(detected_taxa)[1]>1) {
    plist <- list(p_1)
    for(q in 2:dim(detected_taxa)[1]) {
      plist[[q]] <- get(paste0("p_",q))
    }
  } else if (dim(detected_taxa)[1]>0) {
    plist <- list(p_1)
  }
  write.csv(otu5, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".csv",sep=""))
  
  
  ##plot them out
  # tiff(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".tiff",sep=""), units = "px", pointsize = 12, compression = "none", res = 100)
  # do.call(grid.arrange, c(plist))
  # dev.off()
  
  jpeg(h=800*1.25, w=1600*1.25, paste("ancom_",var1,"_",file_path,"_",taxonomic_level,".jpg",sep=""), units = "px", quality = 0.9)
  do.call(grid.arrange, c(plist))
  dev.off()
}
