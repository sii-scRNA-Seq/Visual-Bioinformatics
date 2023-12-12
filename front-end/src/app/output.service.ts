import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { BehaviorSubject, Observable, firstValueFrom } from 'rxjs';
import { Output } from './output';
import { Block } from './block.service';

@Injectable({
  providedIn: 'root',
})

export class OutputService {
  private readonly outputs$: BehaviorSubject<Output[]> = new BehaviorSubject<Output[]> ([]);
  readonly outputs: Observable<Output[]> = this.outputs$.asObservable();

  constructor(private http: HttpClient) { }

  resetOutputs(): void {
    this.outputs$.next([]);
  }
  
  async executeBlock(block: Block): Promise<void> {
    let url = 'http://127.0.0.1:5000/' + block.blockId + '/';
    for (let i=0; i < block.parameters.length; i++) {
      url = url + block.parameters[i].value + '/';
    }
    let outputResponse: Output = {text: '', other: ''};
    try {
      outputResponse = await firstValueFrom(this.http.get<Output>(url));
    } catch(e) {
      outputResponse = {
        text: 'Our servers are not available to complete this request. Please try again. If the issue persists, ...',
        other: ''
      };
    } finally {
      const outputs = this.outputs$.getValue();
      outputs.push(outputResponse);
      this.outputs$.next(outputs);
    }
  }
}
