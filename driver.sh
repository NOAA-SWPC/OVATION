while true
do
  echo "Model and plot generation starting at: $(date)"
  cd /source
  python ovation_model.py && cd /source_plots && python ovation_plot.py
  echo "Model and plot generation completed at: $(date)"
  find /Output/* -mtime +5 -type f -exec rm {} \;
  sleep 600
done
