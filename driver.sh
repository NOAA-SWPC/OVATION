while true
do
  echo "Model and plot generation starting at: $(date)"
  cd /source
  python ovation_model.py && cd /source_plots && python ovation_plot.py
  echo "Model and plot generation completed at: $(date)"
  sleep 600
done
