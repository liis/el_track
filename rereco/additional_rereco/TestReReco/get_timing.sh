cmsRun makeTrackValTree_reTracking.py 1>part1.log 2>part2.log
grep TimeModule\> part2.log > TimingInfo.txt
g++ timing.cpp
./a.out 
