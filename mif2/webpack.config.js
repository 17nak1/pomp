var path = require('path');

module.exports = {
  entry: './src/index.js',
  output: {
    path: path.resolve(__dirname, 'www'),
    filename: 'js/bundle.js'
  },
  node: {
    fs: 'empty',
    path: 'empty'
  },
  devServer: {
    contentBase: path.join(__dirname, 'www'),
    publicPath: '/js/',
    watchContentBase: true,
    compress: true,
    port: 8080
  }
}; 