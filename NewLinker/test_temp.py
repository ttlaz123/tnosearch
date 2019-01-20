from LinkerLib import setTrackId as std
from LinkerLib import writeDetToOrb as std2
from LinkerLib import Detection

num = "chunk1234"
savename = "firstpart+secondpart+iknowcounting+alriteletsgetserious+chunk1234+SNFAKE_SEASON240_ML10.bitch"
print (std(num, savename))

det11 = Detection(12.3, 23.4, 82532.9283, 3722112.23, 12345, 12345)
det12 = Detection(23.4, 23.322, 2323.323123, 23123123.2232, 12346, 23456)
det21 = Detection(123.4, 21.1322, 1222332.2, 21231232.231, 23456, 23432)
tracks = {312819388268398: [det11, det12], 768375293862983: [det21]}
std2 (tracks, "stdOut.dat")
