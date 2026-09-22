import { Component, ErrorInfo, ReactNode } from 'react';

interface Props {
  children?: ReactNode;
  fallback?: ReactNode;
}

interface State {
  hasError: boolean;
  error: Error | null;
}

export default class ErrorBoundary extends Component<Props, State> {
  public state: State = {
    hasError: false,
    error: null
  };

  public static getDerivedStateFromError(error: Error): State {
    return { hasError: true, error };
  }

  public componentDidCatch(error: Error, errorInfo: ErrorInfo) {
    console.error("ErrorBoundary caught an error:", error, errorInfo);
  }

  public render() {
    if (this.state.hasError) {
      return this.props.fallback || (
        <div style={{
          padding: '1rem',
          border: '1px solid #e8542a',
          borderRadius: '6px',
          backgroundColor: '#1c1310',
          color: '#e8542a',
          margin: '0.75rem 0',
          fontSize: '0.9rem',
          fontFamily: 'monospace'
        }}>
          <strong style={{ display: 'block', marginBottom: '0.4rem' }}>⚠️ Error de Renderizado:</strong>
          <span style={{ fontSize: '0.82rem', color: '#fcdcd4' }}>{this.state.error?.message || 'Error desconocido'}</span>
        </div>
      );
    }

    return this.props.children;
  }
}
