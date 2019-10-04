module.exports = {
    verbose: true,
    coverageThreshold: {
        global: {
            statements: 90,
            branches: 90,
            functions: 90,
            lines: 90
        },
        file: {
            statements: 90,
            branches: 90,
            functions: 90,
            lines: 90
        }
    },
    testMatch: [
        '**/test/**/*.spec.js'
    ]
};