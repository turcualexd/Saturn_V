function angle = pitch(t)

t_discr = [0.30,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,50.0,51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,60.0,61.0,62.0,63.0,64.0,65.0,66.0,66.50,67.0,68.0,69.0,70.0,71.0,72.0,73.0,74.0,75.0,76.0,77.0,78.0,79.0,80.0,81.0,82.0,83.0,83.70,84.0,85.0,86.0,87.0,88.0,89.0,90.0,91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100.0,101.0,102.0,103.0,104.0,105.0,106.0,107.0,108.0,109.0,110.0,111.0,112.0,113.0,114.0,115.0,116.0,117.0,118.0,119.0,120.0,121.0,122.0,123.0,124.0,125.0,126.0,127.0,128.0,129.0,130.0,131.0,132.0,133.0,134.0,135.0,135.20,136.0,137.0,138.0,139.0,140.0,141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0,151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,160.0,161.63];
angle_discr = [0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.73,1.46,2.18,2.91,3.64,4.37,5.10,5.82,6.55,7.28,8.01,8.74,9.46,10.19,10.92,11.65,12.38,13.10,13.83,14.56,15.29,16.02,16.74,17.47,18.20,18.93,19.66,20.38,21.11,21.84,22.57,23.30,24.02,24.75,25.48,26.21,26.57,26.94,27.66,28.39,29.12,29.85,30.58,31.30,32.03,32.76,33.49,34.22,34.94,35.67,36.40,36.87,37.34,37.81,38.14,38.28,38.75,39.22,39.69,40.16,40.63,41.10,41.57,42.04,42.51,42.97,43.44,43.91,44.38,44.85,45.32,45.79,46.26,46.73,47.20,47.67,48.14,48.61,49.08,49.55,50.02,50.49,50.96,51.43,51.90,52.37,52.84,53.31,53.78,54.25,54.72,55.19,55.66,56.12,56.59,57.06,57.53,58.00,58.47,58.94,59.41,59.88,60.35,60.82,61.29,61.76,62.23,62.29,62.53,62.82,63.12,63.42,63.72,64.01,64.31,64.61,64.90,65.20,65.50,65.79,66.09,66.39,66.69,66.98,67.28,67.58,67.87,68.17,68.47,68.76,69.06,69.36,69.66,70.14];
angle_discr = deg2rad(angle_discr);

angle = abs(interp1(t_discr, angle_discr, t, 'spline'));