var path = require('path');

module.exports = {
  mode: 'production',
  entry: './src/worker-index.js',
  output: {
    path: path.resolve(__dirname, "../mif2/www", "js"),
    filename: 'worker-bundle-pfilter.js'
  },
  output: {
    path: path.resolve(__dirname, "www", "js"),
    filename: 'worker-bundle.js'
  },
  node: {
    fs: 'empty',
    path: 'empty'
  }
};