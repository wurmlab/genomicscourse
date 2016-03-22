var gulp = require('gulp');
var markdown = require('gulp-markdown');
var reveal = require('gulp-reveal');
var server = require('gulp-server-livereload');

gulp.task('default', function () {
  gulp.src('index.md')
    .pipe(markdown())
    .pipe(reveal())
    .pipe(gulp.dest('.'));
});

gulp.task('watch', function () {
  gulp.watch('*.md', ['default']);
});

gulp.task('serve', ['default', 'watch'], function() {
  gulp.src('.')
    .pipe(server({
      livereload: true,
      directoryListing: false,
      open: true
    }));
});
