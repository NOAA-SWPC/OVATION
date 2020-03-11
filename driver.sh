while true
do
  cd /source
  python ovation_model.py && cd /source_plots && python ovation_plot.py
  echo "Model and plot generation complete"
  sleep 600
done